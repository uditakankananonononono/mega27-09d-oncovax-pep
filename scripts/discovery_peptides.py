"""Discovery layer: aggregate missense mutations from bundle checkpoints,
map to UniProt canonical sequence (QC: WT residue at HGVSp position must
match), build the mutated +/-7 window and all 8-15mer substrings containing
the mutated residue, with WT counterparts. Addendum funnel steps 1-2.
python3 scripts/discovery_peptides.py
"""
import gzip, glob, json, re

IDX = "results/discovery/human_proteome_index.json.gz"
OUT = "results/discovery/neopeptides.json.gz"
AA = set("ACDEFGHIKLMNPQRSTVWY")

def windows(seq, pos0):
    """All (start, length) 8-15mer substrings of seq containing pos0."""
    out = []
    for L in range(8, 16):
        lo = max(0, pos0 - L + 1)
        hi = min(pos0, len(seq) - L)
        for s in range(lo, hi + 1):
            out.append(seq[s:s + L])
    return sorted(set(out))

def main():
    idx = json.load(gzip.open(IDX, "rt"))
    muts, seen = [], set()
    for ck in sorted(glob.glob("results/discovery/maf_bundles/bundle_*.json")):
        b = json.load(open(ck))
        for m in b["mutations"]:
            key = (m["sample"], m["gene"], m["hgvsp"])
            if key not in seen:
                seen.add(key)
                muts.append(m)
    print(f"{len(muts)} unique (sample,gene,hgvsp) missense mutations", flush=True)
    no_gene, qc_fail, ok, recs = 0, 0, 0, []
    for m in muts:
        seq = idx.get(m["gene"])
        if seq is None:
            no_gene += 1
            continue
        p0 = m["pos"] - 1
        if p0 >= len(seq) or seq[p0] != m["wt_aa"]:
            qc_fail += 1
            continue
        if m["mut_aa"] not in AA:
            qc_fail += 1
            continue
        mutseq = seq[:p0] + m["mut_aa"] + seq[p0 + 1:]
        peps = windows(mutseq, p0)
        wtpeps = windows(seq, p0)
        ok += 1
        recs.append({"sample": m["sample"], "gene": m["gene"], "hgvsp": m["hgvsp"],
                     "pos": m["pos"], "wt_aa": m["wt_aa"], "mut_aa": m["mut_aa"],
                     "peptides": peps, "wt_peptides": wtpeps})
    print(f"mapped {ok}; gene-not-in-proteome {no_gene}; QC-fail (WT mismatch) {qc_fail}", flush=True)
    with gzip.open(OUT, "wt") as fh:
        json.dump({"qc": {"n_mutations": len(muts), "mapped": ok,
                          "gene_not_in_proteome": no_gene, "wt_mismatch": qc_fail},
                   "records": recs}, fh)
    print("PEPTIDES DONE", flush=True)

if __name__ == "__main__":
    main()
