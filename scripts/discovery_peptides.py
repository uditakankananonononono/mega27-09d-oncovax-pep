"""Discovery layer (STREAMING, 2GB-box safe): aggregate missense mutations
from bundle checkpoints, map to UniProt canonical sequence (WT-residue QC),
emit mutated +/-7 windows (all 8-15mers containing the mutated residue) +
WT counterparts as JSONL.gz - one mutation per line, never all in memory.
Addendum funnel steps 1-2. python3 scripts/discovery_peptides.py
"""
import gzip, glob, json

IDX = "results/discovery/human_proteome_index.json.gz"
OUT = "results/discovery/neopeptides.jsonl.gz"
AA = set("ACDEFGHIKLMNPQRSTVWY")

def windows(seq, pos0):
    out = []
    for L in range(8, 16):
        lo = max(0, pos0 - L + 1)
        hi = min(pos0, len(seq) - L)
        for s in range(lo, hi + 1):
            out.append(seq[s:s + L])
    return sorted(set(out))

def main():
    idx = json.load(gzip.open(IDX, "rt"))
    seen = set()
    n_mut = n_map = no_gene = qc_fail = 0
    with gzip.open(OUT, "wt") as out:
        for ck in sorted(glob.glob("results/discovery/maf_bundles/bundle_*.json")):
            b = json.load(open(ck))
            for m in b["mutations"]:
                key = (m["sample"], m["gene"], m["hgvsp"])
                if key in seen:
                    continue
                seen.add(key)
                n_mut += 1
                seq = idx.get(m["gene"])
                if seq is None:
                    no_gene += 1
                    continue
                p0 = m["pos"] - 1
                if p0 >= len(seq) or seq[p0] != m["wt_aa"] or m["mut_aa"] not in AA:
                    qc_fail += 1
                    continue
                mutseq = seq[:p0] + m["mut_aa"] + seq[p0 + 1:]
                rec = {"sample": m["sample"], "gene": m["gene"], "hgvsp": m["hgvsp"],
                       "pos": m["pos"], "wt_aa": m["wt_aa"], "mut_aa": m["mut_aa"],
                       "peptides": windows(mutseq, p0), "wt_peptides": windows(seq, p0)}
                out.write(json.dumps(rec) + "\n")
                n_map += 1
            del b
            print(f"{ck}: cumulative {n_map} mapped", flush=True)
    qc = {"n_mutations_unique": n_mut, "mapped": n_map,
          "gene_not_in_proteome": no_gene, "wt_mismatch_or_bad_aa": qc_fail}
    json.dump(qc, open("results/discovery/neopeptides_qc.json", "w"), indent=1)
    print("PEPTIDES DONE", qc, flush=True)

if __name__ == "__main__":
    main()
