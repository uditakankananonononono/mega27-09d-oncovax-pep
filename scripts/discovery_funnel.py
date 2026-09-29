"""Discovery layer: locked funnel (addendum docs/ADDENDUM_DISCOVERY_EXEC_2026-09-28.md).
STREAMING rewrite for the 2GB box (JSONL in/out, external sort, sharded Aho,
disk-chunked mhcflurry cache). Phases: expr -> self -> prescr -> score -> thresh.
python3 scripts/discovery_funnel.py [phase|all]
"""
import csv, gzip, json, math, os, pickle, subprocess, sys
from collections import Counter

PANEL = ["HLA-A*02:01","HLA-A*02:03","HLA-A*03:01","HLA-A*11:01",
         "HLA-A*31:01","HLA-A*68:02","HLA-B*07:02","HLA-B*15:01"]
D = "results/discovery"
KD = dict(A=1.8,C=2.5,D=-3.5,E=-3.5,F=2.8,G=-0.4,H=-3.2,I=4.5,K=-3.9,L=3.8,
          M=1.9,N=-3.5,P=-1.6,Q=-3.5,R=-4.5,S=-0.8,T=-0.7,V=4.2,W=-0.9,Y=-1.3)
CH = dict(D=-1,E=-1,K=1,R=1,H=0.1)
ARO = set('FWY'); AAS = 'ACDEFGHIKLMNPQRSTVWY'
SHARD = 200000
CHUNK = 20000

def phase_expr():
    tpm_cache = {}
    def tpm(sample, gene):
        if sample not in tpm_cache:
            p = f"{D}/expression/{sample}.json"
            tpm_cache[sample] = json.load(open(p))["tpm"] if os.path.exists(p) else None
            if len(tpm_cache) > 50:      # bound cache
                tpm_cache.pop(next(iter(tpm_cache)))
        t = tpm_cache[sample]
        return None if t is None else t.get(gene)
    n_in = n_surv = no_expr = low = 0
    with gzip.open(f"{D}/neopeptides.jsonl.gz", "rt") as fin, \
         gzip.open(f"{D}/funnel_expr.jsonl.gz", "wt") as fout:
        for line in fin:
            r = json.loads(line)
            n_in += 1
            t = tpm(r["sample"], r["gene"])
            if t is None:
                no_expr += 1
                continue
            if t < 1.0:
                low += 1
                continue
            r["tpm"] = t
            fout.write(json.dumps(r) + "\n")
            n_surv += 1
    json.dump({"n_in": n_in, "n_surv": n_surv, "no_expression_data": no_expr,
               "below_1tpm": low}, open(f"{D}/funnel_expr_qc.json", "w"), indent=1)
    print(f"expr: {n_surv}/{n_in} survive (no-expr {no_expr}, <1TPM {low})", flush=True)

def _unique_peptides(infile, outfile):
    """External-sort dedup of every peptide instance; RAM-lean."""
    with gzip.open(infile, "rt") as fin, open(f"{D}/_peps.txt", "w") as tmp:
        for line in fin:
            for p in json.loads(line)["peptides"]:
                tmp.write(p + "\n")
    subprocess.run(f"sort -u {D}/_peps.txt | gzip > {outfile}", shell=True, check=True)
    os.remove(f"{D}/_peps.txt")

def phase_self():
    import ahocorasick
    _unique_peptides(f"{D}/funnel_expr.jsonl.gz", f"{D}/_uniq_peps.txt.gz")
    idx = json.load(gzip.open(f"{D}/human_proteome_index.json.gz", "rt"))
    seqs = list(idx.values())
    del idx
    self_hits = set()
    shard, n_tot = [], 0
    def scan(shard):
        A = ahocorasick.Automaton()
        for p in shard:
            A.add_word(p, p)
        A.make_automaton()
        hits = set()
        for seq in seqs:
            for _, p in A.iter(seq):
                hits.add(p)
        return hits
    with gzip.open(f"{D}/_uniq_peps.txt.gz", "rt") as fh:
        for line in fh:
            shard.append(line.strip())
            n_tot += 1
            if len(shard) >= SHARD:
                self_hits |= scan(shard)
                print(f"self: scanned {n_tot} peptides, {len(self_hits)} self-hits so far", flush=True)
                shard = []
    if shard:
        self_hits |= scan(shard)
    del seqs
    n_in = n_surv = 0
    with gzip.open(f"{D}/funnel_expr.jsonl.gz", "rt") as fin, \
         gzip.open(f"{D}/funnel_self.jsonl.gz", "wt") as fout:
        for line in fin:
            r = json.loads(line)
            n_in += 1
            keep = [p for p in r["peptides"] if p not in self_hits]
            if keep:
                r["peptides"] = keep
                fout.write(json.dumps(r) + "\n")
                n_surv += 1
    json.dump({"n_unique_peptides": n_tot, "n_self_excluded": len(self_hits),
               "n_records_in": n_in, "n_records_surv": n_surv},
              open(f"{D}/funnel_self_qc.json", "w"), indent=1)
    print(f"self: {len(self_hits)}/{n_tot} unique peptides self-excluded; {n_surv}/{n_in} records survive", flush=True)

def _mhcflurry(pred, pairs):
    res = pred.predict_to_dataframe(peptides=[p for p, a in pairs],
                                    alleles=[a for p, a in pairs], throw=False)
    out = []
    for r in res.itertuples():
        pct = getattr(r, "prediction_percentile", None)
        out.append([r.peptide, r.allele, float(r.prediction),
                    None if pct is None else float(pct)])
    return out

def phase_prescr():
    # STREAMING (2GB box): never materialize the full peptide/pair lists.
    # Chunk boundaries identical to the old pep-major pairs list: each chunk =
    # CHUNK/len(PANEL) consecutive peptides from _uniq_self.txt.gz x all alleles.
    from mhcflurry import Class1AffinityPredictor
    os.makedirs(f"{D}/mhc_chunks", exist_ok=True)
    _unique_peptides(f"{D}/funnel_self.jsonl.gz", f"{D}/_uniq_self.txt.gz")
    n_peps = sum(1 for _ in gzip.open(f"{D}/_uniq_self.txt.gz", "rt"))
    n_pairs = n_peps * len(PANEL)
    print(f"prescr: {n_peps} unique peptides x {len(PANEL)} = {n_pairs} pairs", flush=True)
    pred = Class1AffinityPredictor.load()
    done_chunks = {int(f.split('_')[1].split('.')[0]) for f in os.listdir(f"{D}/mhc_chunks")}
    per_chunk = CHUNK // len(PANEL)
    buf, ci = [], 0
    done_this_session = [0]
    REEXEC_EVERY = 60  # re-exec fresh interpreter: RSS growth slows TF ~8x by chunk ~30
    def flush(buf, ci):
        pairs = [(p, a) for p in buf for a in PANEL]
        out = _mhcflurry(pred, pairs)
        json.dump(out, gzip.open(f"{D}/mhc_chunks/chunk_{ci//CHUNK:05d}.json.gz", "wt"))
        done_this_session[0] += 1
        if (ci // CHUNK) % 25 == 0:
            print(f"  mhcflurry chunk {ci//CHUNK} ({ci}/{n_pairs})", flush=True)
        if done_this_session[0] >= REEXEC_EVERY:
            print("prescr: re-exec for fresh interpreter", flush=True)
            os.execv(sys.executable, [sys.executable] + sys.argv)
    box_end = int(os.environ.get("PRESCR_BOX_END_CHUNK", "999999"))
    for line in gzip.open(f"{D}/_uniq_self.txt.gz", "rt"):
        p = line.strip()
        if not p: continue
        if ci // CHUNK >= box_end:
            print(f"prescr: stopped before shard-assigned chunk {box_end}", flush=True)
            return
        if ci // CHUNK in done_chunks:
            ci += len(PANEL)
            continue
        buf.append(p)
        if len(buf) == per_chunk:
            flush(buf, ci)
            buf = []
        ci += len(PANEL)
    if buf:
        # final partial chunk (skip only if its index is done)
        if ci // CHUNK not in done_chunks:
            flush(buf, ci)
    # streaming pass: peptides with any pct<=2
    keep_peps = set()
    for f in sorted(os.listdir(f"{D}/mhc_chunks")):
        for p, a, aff, pct in json.load(gzip.open(f"{D}/mhc_chunks/{f}", "rt")):
            if pct is not None and pct <= 2.0:
                keep_peps.add(p)
    n_in = n_surv = 0
    with gzip.open(f"{D}/funnel_self.jsonl.gz", "rt") as fin, \
         gzip.open(f"{D}/funnel_prescr.jsonl.gz", "wt") as fout:
        for line in fin:
            r = json.loads(line)
            n_in += 1
            kp = [p for p in r["peptides"] if p in keep_peps]
            if kp:
                r["peptides"] = kp
                fout.write(json.dumps(r) + "\n")
                n_surv += 1
    json.dump({"n_pairs": n_pairs, "n_unique_peptides": n_peps,
               "n_peptides_pass_pct2": len(keep_peps),
               "n_records_in": n_in, "n_records_surv": n_surv},
              open(f"{D}/funnel_prescr_qc.json", "w"), indent=1)
    print(f"prescr: {len(keep_peps)} peptides pass pct<=2; {n_surv}/{n_in} records survive", flush=True)

def feats(pep, aff_nm):
    import numpy as np
    s = pep; n = len(s)
    c = Counter(s)
    comp = [c.get(a, 0) / n for a in AAS]
    kd = [KD[x] for x in s]
    return comp + [n, float(np.mean(kd)), float(np.max(kd) - np.min(kd)),
                   sum(CH.get(x, 0) for x in s),
                   sum(1 for x in s if x in ARO) / n,
                   math.log(max(aff_nm, 1e-3))]

def phase_score():
    import numpy as np
    rows = list(csv.DictReader(open("data/processed/splits_v1.csv")))
    tr = [r for r in rows if r["split"] == "train"]
    aff_train = {}
    cache = "/tmp/09d_aff_cache.pkl"
    if os.path.exists(cache):
        aff_train = pickle.load(open(cache, "rb"))
    miss = [(r["peptide"], r["allele"]) for r in tr if (r["peptide"], r["allele"]) not in aff_train]
    if miss:
        from mhcflurry import Class1AffinityPredictor
        pred = Class1AffinityPredictor.load()
        for p, a, aff, pct in _mhcflurry(pred, miss):
            aff_train[(p, a)] = aff
        pickle.dump(aff_train, open(cache, "wb"))
    Xtr = np.array([feats(r["peptide"], aff_train.get((r["peptide"], r["allele"]), 500.0)) for r in tr])
    ytr = np.array([int(r["label"]) for r in tr])
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import HistGradientBoostingClassifier
    lr = LogisticRegression(max_iter=2000, C=1.0).fit(Xtr, ytr)
    hgb = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.08,
                                         max_depth=None, random_state=2709).fit(Xtr, ytr)
    # kept pairs from prescreen records
    keep = set()
    with gzip.open(f"{D}/funnel_prescr.jsonl.gz", "rt") as fin:
        for line in fin:
            for p in json.loads(line)["peptides"]:
                for a in PANEL:
                    keep.add((p, a))
    # affinities for kept pairs from chunks (only kept pairs held in RAM)
    aff = {}
    for f in sorted(os.listdir(f"{D}/mhc_chunks")):
        for p, a, af, pct in json.load(gzip.open(f"{D}/mhc_chunks/{f}", "rt")):
            if (p, a) in keep and pct is not None and pct <= 2.0:
                aff[(p, a)] = af
    pairs = sorted(aff)
    print(f"score: {len(pairs)} prescreened pairs", flush=True)
    X = np.array([feats(p, aff[(p, a)]) for p, a in pairs])
    with gzip.open(f"{D}/funnel_scores.jsonl.gz", "wt") as out:
        for name, clf in [("hgb", hgb), ("logreg", lr)]:
            s = clf.predict_proba(X)[:, 1]
            for (p, a), v in zip(pairs, s):
                out.write(json.dumps({"model": name, "peptide": p, "allele": a,
                                      "score": float(v), "aff_nm": aff[(p, a)]}) + "\n")
    print("score: done", flush=True)

def phase_thresh():
    import numpy as np
    vals = []
    with gzip.open(f"{D}/funnel_scores.jsonl.gz", "rt") as fh:
        for line in fh:
            r = json.loads(line)
            if r["model"] == "hgb":
                vals.append(r["score"])
    vals = np.array(sorted(vals))
    thr = float(np.quantile(vals, 0.99))
    named = {}
    with gzip.open(f"{D}/funnel_scores.jsonl.gz", "rt") as fh:
        for line in fh:
            r = json.loads(line)
            if r["model"] == "hgb" and r["score"] >= thr:
                named[f"{r['peptide']}|{r['allele']}"] = r["score"]
    json.dump({"n_scored_pairs": len(vals), "top1pct_threshold": thr,
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
