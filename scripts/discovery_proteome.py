"""Discovery layer: download UniProt reviewed human proteome (UP000005640)
and build a gene -> canonical-sequence index. Serves both the +/-7 window
context for neopeptide generation and the self-proteome filter
(addendum docs/ADDENDUM_DISCOVERY_EXEC_2026-09-28.md).
One file, hashed; index committed gzipped. python3 scripts/discovery_proteome.py
"""
import gzip, hashlib, json, os, urllib.request

URL = ("https://rest.uniprot.org/uniprotkb/stream?query="
       "%28proteome%3AUP000005640%29%20AND%20%28reviewed%3Atrue%29&format=fasta&compressed=true")
RAW = "results/discovery/up000005640_reviewed.fasta.gz"
IDX = "results/discovery/human_proteome_index.json.gz"

def fetch():
    os.makedirs("results/discovery", exist_ok=True)
    if not os.path.exists(RAW):
        print("downloading proteome...", flush=True)
        urllib.request.urlretrieve(URL, RAW + ".tmp")
        os.rename(RAW + ".tmp", RAW)
    h = hashlib.sha256(open(RAW, "rb").read()).hexdigest()
    print("sha256:", h, flush=True)
    return h

def build():
    idx, amb = {}, []
    name, seq, gn = None, [], None
    def flush():
        if name is None: return
        s = "".join(seq)
        if gn:
            if gn in idx and idx[gn] != s:
                amb.append(gn)   # duplicate gene name, different sequence: keep first, disclose
            else:
                idx.setdefault(gn, s)
    with gzip.open(RAW, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                gn = None
                for tok in line.split():
                    if tok.startswith("GN="):
                        gn = tok[3:]
                        break
                seq = []
            else:
                seq.append(line.strip())
    flush()
    with gzip.open(IDX, "wt") as fh:
        json.dump(idx, fh)
    print(f"indexed {len(idx)} genes, {len(amb)} ambiguous (first kept, disclosed)", flush=True)
    return idx, amb

if __name__ == "__main__":
    sha = fetch()
    idx, amb = build()
    meta = {"source_url": URL, "sha256_fasta_gz": sha, "n_genes": len(idx),
            "ambiguous_gene_names_first_kept": amb[:50], "n_ambiguous": len(amb),
            "downloaded": "2026-09-28"}
    json.dump(meta, open("results/discovery/proteome_meta.json", "w"), indent=1)
    print("PROTEOME DONE", flush=True)
