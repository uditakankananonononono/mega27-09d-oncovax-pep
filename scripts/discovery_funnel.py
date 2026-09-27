"""Discovery layer: locked funnel (addendum docs/ADDENDUM_DISCOVERY_EXEC_2026-09-28.md).
Phases, resume-safe, RAM-staged (2GB box):
  expr   - expression pre-filter (carrier-sample TPM >= 1)
  self   - self-proteome exact-match filter (Aho-Corasick over proteome index)
  prescr - MHCflurry-2.2.1 affinity percentile <= 2.0 vs >=1 panel allele
  score  - v1 retrain (locked hyperparams, train split only) + score survivors
  thresh - top-1% hgb threshold locked in JSON before named selection
python3 scripts/discovery_funnel.py [phase]
"""
import csv, glob, gzip, json, math, os, pickle, sys
from collections import Counter, defaultdict

PANEL = ["HLA-A*02:01","HLA-A*02:03","HLA-A*03:01","HLA-A*11:01",
         "HLA-A*31:01","HLA-A*68:02","HLA-B*07:02","HLA-B*15:01"]
D = "results/discovery"
KD = dict(A=1.8,C=2.5,D=-3.5,E=-3.5,F=2.8,G=-0.4,H=-3.2,I=4.5,K=-3.9,L=3.8,
          M=1.9,N=-3.5,P=-1.6,Q=-3.5,R=-4.5,S=-0.8,T=-0.7,V=4.2,W=-0.9,Y=-1.3)
CH = dict(D=-1,E=-1,K=1,R=1,H=0.1)
ARO = set('FWY'); AAS = 'ACDEFGHIKLMNPQRSTVWY'

def phase_expr():
    pep = json.load(gzip.open(f"{D}/neopeptides.json.gz", "rt"))
    tpm_cache = {}
    def tpm(sample, gene):
        if sample not in tpm_cache:
            p = f"{D}/expression/{sample}.json"
            tpm_cache[sample] = json.load(open(p))["tpm"] if os.path.exists(p) else None
        t = tpm_cache[sample]
        return None if t is None else t.get(gene)
    surv, no_expr, low_expr = [], 0, 0
    for r in pep["records"]:
        t = tpm(r["sample"], r["gene"])
        if t is None:
            no_expr += 1
            continue
        if t >= 1.0:
            r = dict(r); r["tpm"] = t
            surv.append(r)
        else:
            low_expr += 1
    json.dump({"n_in": len(pep["records"]), "n_surv": len(surv),
               "no_expression_data": no_expr, "below_1tpm": low_expr,
               "records": surv}, gzip.open(f"{D}/funnel_expr.json.gz", "wt"))
    print(f"expr: {len(surv)}/{len(pep['records'])} survive (no-expr {no_expr}, <1TPM {low_expr})", flush=True)

def phase_self():
    import ahocorasick
    recs = json.load(gzip.open(f"{D}/funnel_expr.json.gz", "rt"))["records"]
    peps = sorted({p for r in recs for p in r["peptides"]})
    A = ahocorasick.Automaton()
    for i, p in enumerate(peps):
        A.add_word(p, i)
    A.make_automaton()
    self_hits = set()
    idx = json.load(gzip.open(f"{D}/human_proteome_index.json.gz", "rt"))
    for gene, seq in idx.items():
        for _, i in A.iter(seq):
            self_hits.add(peps[i])
    print(f"self: {len(self_hits)}/{len(peps)} peptides exact-match human proteome -> excluded", flush=True)
    out = []
    for r in recs:
        keep = [p for p in r["peptides"] if p not in self_hits]
        if keep:
            r = dict(r); r["peptides"] = keep
            out.append(r)
    json.dump({"n_peptides_in": len(peps), "n_self_excluded": len(self_hits),
               "n_records_surv": len(out), "records": out},
              gzip.open(f"{D}/funnel_self.json.gz", "wt"))
    print(f"self: {len(out)} records survive", flush=True)

def _mhcflurry(pairs):
    from mhcflurry import Class1AffinityPredictor
    pred = Class1AffinityPredictor.load()
    peps = sorted({p for p, a in pairs})
    res = pred.predict_to_dataframe(peptides=[p for p, a in pairs],
                                    alleles=[a for p, a in pairs], throw=False)
    out = {}
    for r in res.itertuples():
        pct = getattr(r, "prediction_percentile", None)
        out[(r.peptide, r.allele)] = {"aff": float(r.prediction),
                                      "pct": None if pct is None else float(pct)}
    return out

def phase_prescr():
    recs = json.load(gzip.open(f"{D}/funnel_self.json.gz", "rt"))["records"]
    peps = sorted({p for r in recs for p in r["peptides"]})
    pairs = [(p, a) for p in peps for a in PANEL]
    print(f"prescr: {len(peps)} peptides x {len(PANEL)} alleles = {len(pairs)} pairs", flush=True)
    cache = f"{D}/mhcflurry_cache.pkl"
    aff = pickle.load(open(cache, "rb")) if os.path.exists(cache) else {}
    todo = [pr for pr in pairs if pr not in aff]
    for i in range(0, len(todo), 20000):
        chunk = todo[i:i + 20000]
        aff.update(_mhcflurry(chunk))
        pickle.dump(aff, open(cache, "wb"))
        print(f"  mhcflurry {min(i+20000,len(todo))}/{len(todo)}", flush=True)
    keep_pairs = {pr for pr, v in aff.items() if v["pct"] is not None and v["pct"] <= 2.0}
    keep_peps = {p for p, a in keep_pairs}
    out = []
    for r in recs:
        kp = [p for p in r["peptides"] if p in keep_peps]
        if kp:
            r = dict(r); r["peptides"] = kp
            r["best_pct"] = {p: min(aff[(p, a)]["pct"] for a in PANEL if (p, a) in aff and aff[(p,a)]["pct"] is not None) for p in kp}
            out.append(r)
    json.dump({"n_pairs": len(pairs), "n_pairs_pass_pct2": len(keep_pairs),
               "n_records_surv": len(out), "records": out},
              gzip.open(f"{D}/funnel_prescr.json.gz", "wt"))
    print(f"prescr: {len(keep_pairs)} pairs pass pct<=2; {len(out)} records survive", flush=True)

def feats(pep, allele, aff):
    import numpy as np
    s = pep; n = len(s)
    c = Counter(s)
    comp = [c.get(a, 0) / n for a in AAS]
    kd = [KD[x] for x in s]
    a = aff.get((s, allele), {}).get("aff", 500.0)
    return comp + [n, float(np.mean(kd)), float(np.max(kd) - np.min(kd)),
                   sum(CH.get(x, 0) for x in s),
                   sum(1 for x in s if x in ARO) / n,
                   math.log(max(a, 1e-3))]

def phase_score():
    import numpy as np, pickle as pk
    aff = pk.load(open(f"{D}/mhcflurry_cache.pkl", "rb"))
    # retrain v1 on frozen train split, locked hyperparams (train_model_v1.py)
    rows = list(csv.DictReader(open("data/processed/splits_v1.csv")))
    tr = [r for r in rows if r["split"] == "train"]
    train_pairs = sorted({(r["peptide"], r["allele"]) for r in tr})
    missing = [pr for pr in train_pairs if pr not in aff]
    if missing:
        aff.update(_mhcflurry(missing))
        pk.dump(aff, open(f"{D}/mhcflurry_cache.pkl", "wb"))
    import numpy as np
    Xtr = np.array([feats(r["peptide"], r["allele"], aff) for r in tr])
    ytr = np.array([int(r["label"]) for r in tr])
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import HistGradientBoostingClassifier
    lr = LogisticRegression(max_iter=2000, C=1.0).fit(Xtr, ytr)
    hgb = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.08,
                                         max_depth=None, random_state=2709).fit(Xtr, ytr)
    recs = json.load(gzip.open(f"{D}/funnel_prescr.json.gz", "rt"))["records"]
    pairs = sorted({(p, a) for r in recs for p in r["peptides"] for a in PANEL})
    # only score pairs that passed the prescreen
    passed = set()
    for r in recs:
        for p in r["peptides"]:
            for a in PANEL:
                if aff.get((p, a), {}).get("pct") is not None and aff[(p, a)]["pct"] <= 2.0:
                    passed.add((p, a))
    pairs = sorted(passed)
    print(f"score: {len(pairs)} prescreened pairs", flush=True)
    X = np.array([feats(p, a, aff) for p, a in pairs])
    out = {}
    for name, clf in [("hgb", hgb), ("logreg", lr)]:
        s = clf.predict_proba(X)[:, 1]
        out[name] = {f"{p}|{a}": float(v) for (p, a), v in zip(pairs, s)}
    json.dump({"n_pairs": len(pairs), "scores": out}, gzip.open(f"{D}/funnel_scores.json.gz", "wt"))
    print("score: done", flush=True)

def phase_thresh():
    import numpy as np
    d = json.load(gzip.open(f"{D}/funnel_scores.json.gz", "rt"))
    hgb = d["scores"]["hgb"]
    vals = np.array(sorted(hgb.values()))
    thr = float(np.quantile(vals, 0.99))
    named = {k: v for k, v in hgb.items() if v >= thr}
    json.dump({"n_scored_pairs": len(hgb), "top1pct_threshold": thr,
               "n_named_candidates": len(named),
               "threshold_rule": "hgb >= 99th percentile of all prescreened pair scores (locked before named selection)"},
              open(f"{D}/funnel_threshold.json", "w"), indent=1)
    json.dump(named, open(f"{D}/named_candidates_scores.json", "w"), indent=1)
    print(f"thresh: top-1% = {thr:.5f}, {len(named)} named candidate pairs", flush=True)

if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"
    if phase in ("all", "expr"): phase_expr()
    if phase in ("all", "self"): phase_self()
    if phase in ("all", "prescr"): phase_prescr()
    if phase in ("all", "score"): phase_score()
    if phase in ("all", "thresh"): phase_thresh()
    print("FUNNEL DONE", flush=True)
