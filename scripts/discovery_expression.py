"""Discovery layer: GDC STAR gene counts for mutation-carrying samples.
Maps file -> 16-char tumor sample barcode, downloads only samples with
QC-passing candidate mutations, streams TPM rows for mutated genes,
deletes raw. Addendum funnel step 3 (>= 1 TPM pre-filter applied later,
in the funnel script, from these per-sample values). Resume-safe.
python3 scripts/discovery_expression.py
"""
import gzip, io, json, os, tarfile, time, urllib.request

OUT = "results/discovery/expression"
MAN = "results/discovery/expression_manifest.json"

def post(url, payload, timeout=180):
    req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode()

def build_manifest():
    filt = {"op":"and","content":[
        {"op":"=","content":{"field":"cases.project.project_id","value":"TCGA-SKCM"}},
        {"op":"=","content":{"field":"data_type","value":"Gene Expression Quantification"}},
        {"op":"=","content":{"field":"access","value":"open"}}]}
    fields = "file_id,file_name,md5sum,cases.samples.submitter_id,cases.samples.sample_type"
    out, from_ = [], 0
    while True:
        d = json.loads(post("https://api.gdc.cancer.gov/files",
                            {"filters":filt,"fields":fields,"size":500,"from":from_,"format":"json"}))
        hits = d["data"]["hits"]
        for h in hits:
            for c in h.get("cases", []):
                for s in c.get("samples", []):
                    st = s.get("sample_type", "")
                    if "Tumor" in st:
                        out.append({"file_id": h["file_id"], "file_name": h["file_name"],
                                    "md5sum": h.get("md5sum"),
                                    "sample": s["submitter_id"][:16], "sample_type": st})
        from_ += len(hits)
        if from_ >= d["data"]["pagination"]["total"] or not hits: break
    json.dump({"retrieved": "2026-09-28", "n": len(out), "files": out}, open(MAN, "w"), indent=1)
    return out

def parse_star(b, genes):
    tpm = {}
    text = b.decode("utf-8", "replace")
    hdr = None
    for line in text.split("\n"):
        if line.startswith("#") or not line: continue
        if hdr is None:
            hdr = line.split("\t")
            try:
                ig, it = hdr.index("gene_name"), hdr.index("tpm_unstranded")
            except ValueError:
                return tpm
            continue
        f = line.split("\t")
        if len(f) <= max(ig, it): continue
        if f[ig] in genes:
            try: tpm[f[ig]] = float(f[it])
            except ValueError: pass
    return tpm

def main():
    man = build_manifest() if not os.path.exists(MAN) else json.load(open(MAN))["files"]
    need = {}
    with gzip.open("results/discovery/neopeptides.jsonl.gz", "rt") as fh:
        for line in fh:
            r = json.loads(line)
            need.setdefault(r["sample"], set()).add(r["gene"])
    fmap = {}
    for f in man:
        fmap.setdefault(f["sample"], f)   # first tumor aliquot wins; disclosed
    os.makedirs(OUT, exist_ok=True)
    missing_expr = []
    for i, (sample, genes) in enumerate(sorted(need.items())):
        ck = f"{OUT}/{sample}.json"
        if os.path.exists(ck): continue
        f = fmap.get(sample)
        if f is None:
            missing_expr.append(sample)
            continue
        try:
            req = urllib.request.Request(f"https://api.gdc.cancer.gov/data/{f['file_id']}")
            with urllib.request.urlopen(req, timeout=300) as r:
                b = r.read()
        except Exception as e:
            print(f"{sample}: download failed {e}", flush=True)
            time.sleep(5)
            continue
        if b[:2] == b"\x1f\x8b":
            b = gzip.decompress(b)
        tpm = parse_star(b, genes)
        json.dump({"sample": sample, "file_id": f["file_id"], "file_name": f["file_name"],
                   "md5sum": f["md5sum"], "tpm": tpm}, open(ck, "w"))
        if i % 25 == 0:
            print(f"[{i}/{len(need)}] {sample}: {len(tpm)}/{len(genes)} genes", flush=True)
        time.sleep(0.5)
    json.dump(missing_expr, open("results/discovery/expression_missing.json", "w"))
    print(f"EXPRESSION DONE (missing mapping for {len(missing_expr)} samples, disclosed)", flush=True)

if __name__ == "__main__":
    main()
