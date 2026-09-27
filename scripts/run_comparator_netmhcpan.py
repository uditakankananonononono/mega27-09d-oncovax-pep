#!/usr/bin/env python3
"""Comparator A: NetMHCpan-4.1 EL percentile rank on the frozen test set,
via IEDB nextgen tools API. Addendum docs/ADDENDUM_COMPARATOR_A_EXEC_2026-09-28.md.

Phase 1 (fetch): per (allele, exact length) batch; each peptide its own
FASTA record; peptide_length_range [L,L]; raw API JSON cached in
results/comparator_a_raw/. Resumable - skips existing cache files.
Phase 2 (--score): assemble per-pair scores, symmetric-exclusion report,
metrics via verbatim v1 metric code, write results/comparator_a_netmhcpan.json.
"""
import csv, json, math, os, sys, time, urllib.request

class UnsupportedError(Exception):
    pass

API = "https://api-nextgen-tools.iedb.org/api/v1"
RAW = "results/comparator_a_raw"
SPLITS = "data/processed/splits_v1.csv"

def post_json(url, payload, timeout=120):
    req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())

def get_json(url, timeout=120):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return json.loads(r.read().decode())

def run_batch(allele, L, peps):
    fasta = "".join(f">p{i}\n{p}\n" for i, p in enumerate(peps))
    payload = {"pipeline_id": "", "run_stage_range": [1, 1], "stages": [
        {"stage_number": 1, "tool_group": "mhci",
         "input_sequence_text": fasta,
         "input_parameters": {"alleles": allele,
                              "peptide_length_range": [L, L],
                              "predictors": [{"type": "binding",
                                              "method": "netmhcpan_el"}]}}]}
    sub = post_json(f"{API}/pipeline", payload)
    if "results_uri" not in sub:
        raise UnsupportedError(str(sub.get("errors", sub))[:300])
    uri = sub["results_uri"]
    for _ in range(60):  # up to ~5 min per batch
        time.sleep(5)
        res = get_json(uri)
        if res.get("status") == "done":
            return res
        if res.get("status") not in ("pending", "running", None):
            raise RuntimeError(f"bad status {res.get('status')}")
    raise TimeoutError(f"{allele} L{L} still pending after 5 min")

def fetch_all():
    rows = list(csv.DictReader(open(SPLITS)))
    test = [r for r in rows if r["split"] == "test"]
    alleles = sorted({r["allele"] for r in test})
    if "--shard" in sys.argv:
        i = sys.argv.index("--shard")
        idx, n = int(sys.argv[i + 1]), int(sys.argv[i + 2])
        alleles = [a for j, a in enumerate(alleles) if j % n == idx]
        print(f"shard {idx}/{n}: {len(alleles)} alleles", flush=True)
    by_len = {}
    for r in test:
        by_len.setdefault(len(r["peptide"]), set()).add(r["peptide"])
    os.makedirs(RAW, exist_ok=True)
    total = sum(1 for a in alleles for L in by_len)
    done = 0
    for a in alleles:
        for L, peps in sorted(by_len.items()):
            done += 1
            cache = f"{RAW}/{a.replace('*','')}_L{L}.json"
            if os.path.exists(cache):
                print(f"[{done}/{total}] {a} L{L} cached", flush=True)
                continue
            peps = sorted(peps)
            try:
                res = run_batch(a, L, peps)
            except UnsupportedError as e:
                json.dump({"unsupported": True, "allele": a, "L": L,
                           "api_error": str(e)}, open(cache, "w"))
                print(f"[{done}/{total}] {a} L{L} UNSUPPORTED by API - cached marker", flush=True)
                continue
            except Exception as e:
                ok = False
                for attempt in range(3):
                    try:
                        res = run_batch(a, L, peps)
                        ok = True
                        break
                    except UnsupportedError as e2:
                        res = e2
                        break
                    except Exception as e2:
                        print(f"[{done}/{total}] {a} L{L} attempt {attempt+1} failed: {e2}", flush=True)
                        time.sleep([10, 30, 90][attempt])
                if isinstance(res, UnsupportedError):
                    json.dump({"unsupported": True, "allele": a, "L": L,
                               "api_error": str(res)}, open(cache, "w"))
                    print(f"[{done}/{total}] {a} L{L} UNSUPPORTED by API - cached marker", flush=True)
                elif not ok:
                    print(f"[{done}/{total}] {a} L{L} FAILED 3x - disclosed exclusion", flush=True)
                else:
                    json.dump(res, open(cache, "w"))
                    print(f"[{done}/{total}] {a} L{L} ok ({len(peps)} peps)", flush=True)
                time.sleep(2)
                continue
            json.dump(res, open(cache, "w"))
            print(f"[{done}/{total}] {a} L{L} ok ({len(peps)} peps)", flush=True)
            time.sleep(2)

def score():
    import numpy as np
    from sklearn.metrics import average_precision_score, roc_auc_score

    # --- verbatim v1 metric code (train_model_v1.py), disclosed ---
    def ndcg_at_k(score, gain, k=100):
        order = np.argsort(-score)[:k]
        dcg = sum((2**gain[i] - 1) / math.log2(r + 2) for r, i in enumerate(order))
        ideal = np.argsort(-gain)[:k]
        idcg = sum((2**gain[i] - 1) / math.log2(r + 2) for r, i in enumerate(ideal))
        return dcg / idcg if idcg > 0 else 0.0

    rows = list(csv.DictReader(open(SPLITS)))
    test = [r for r in rows if r["split"] == "test"]
    pct = {}   # (peptide, allele) -> percentile
    missing_batches = []
    for fn in os.listdir(RAW):
        if not fn.endswith(".json"):
            continue
        d = json.load(open(f"{RAW}/{fn}"))
        if d.get("unsupported"):
            continue
        for t in d["data"]["results"]:
            if t["type"] != "peptide_table":
                continue
            cols = [c["name"] for c in t["table_columns"]]
            ip, ia = cols.index("peptide"), cols.index("allele")
            ipct = cols.index("netmhcpan_el_percentile")
            for row in t["table_data"]:
                pct[(row[ip], row[ia])] = float(row[ipct])
    alleles_needed = sorted({r["allele"] for r in test})
    lens_needed = sorted({len(r["peptide"]) for r in test})
    for a in alleles_needed:
        for L in lens_needed:
            if not os.path.exists(f"{RAW}/{a.replace('*','')}_L{L}.json"):
                missing_batches.append(f"{a} L{L}")
    scored, missing_pairs = [], []
    for r in test:
        key = (r["peptide"], r["allele"])
        if key in pct:
            scored.append((r, pct[key]))
        else:
            missing_pairs.append(key)
    print(f"scored {len(scored)}/{len(test)} pairs; missing {len(missing_pairs)}; "
          f"missing batches: {missing_batches}", flush=True)
    transient_missing = [k for k in missing_pairs if len(k[0]) != 15]
    if transient_missing:
        json.dump({"missing_pairs": transient_missing, "missing_batches": missing_batches},
                  open("results/comparator_a_missing.json", "w"), indent=1)
        print("INCOMPLETE (non-L15 gaps) - fix before metrics", flush=True)
        return
    # symmetric exclusion: L15 pairs dropped from BOTH sides (addendum note 1)
    excluded = set(missing_pairs)
    print(f"symmetric exclusion: {len(excluded)} pairs (L15, unsupported by NetMHCpan)", flush=True)
    yte = np.array([int(r["label"]) for r, _ in scored])
    gte = np.array([float(r["positive_fraction"]) for r, _ in scored])
    s = np.array([-p for _, p in scored])  # larger = stronger binder
    k = min(100, len(s))
    order = np.argsort(-s)[:k]
    out = dict(AUPRC=float(average_precision_score(yte, s)),
               AUROC=float(roc_auc_score(yte, s)),
               Recall_at_100=float(yte[order].mean()),
               NDCG_at_100=float(ndcg_at_k(s, gte, k)))
    res = dict(comparator="NetMHCpan-4.1 EL percentile (IEDB nextgen API)",
               addendum="docs/ADDENDUM_COMPARATOR_A_EXEC_2026-09-28.md",
               n_test=len(scored),
               n_excluded_l15=len(excluded),
               excluded_l15=sorted([list(k) for k in excluded]),
               test_pos_rate=float(yte.mean()),
               results={"netmhcpan_el_rank": out},
               model_v1_reference_full_test=json.load(open("results/model_v1_frozen_test.json"))["results"],
               model_v1_reference_subset_8_14={k: v["subset_8_14"] for k, v in
                   json.load(open("results/model_v1_rescore_subset.json"))["results"].items()})
    json.dump(res, open("results/comparator_a_netmhcpan.json", "w"), indent=1)
    print(json.dumps(out, indent=1), flush=True)
    print("SCORE DONE", flush=True)

if __name__ == "__main__":
    if "--score" in sys.argv:
        score()
    else:
        fetch_all()
