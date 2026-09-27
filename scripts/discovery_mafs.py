"""Discovery layer: download TCGA-SKCM masked MAFs per committed manifest,
stream-parse missense SNVs (Variant_Classification == Missense_Mutation,
valid HGVSp_Short p.XnY), aggregate, delete raws (disk). Bundled via GDC
/data POST (ids list -> tar). Resume-safe: bundle checkpoint files.
Addendum docs/ADDENDUM_DISCOVERY_EXEC_2026-09-28.md.
python3 scripts/discovery_mafs.py
"""
import io, json, os, re, tarfile, time, urllib.request

MAN = "results/discovery/skcm_maf_manifest.json"
OUTDIR = "results/discovery/maf_bundles"
BUNDLE = 47          # files per POST bundle (~10 bundles total)
HGVSP = re.compile(r"^p\.([A-Z])(\d+)([A-Z])$")

def post_data(ids, retries=4):
    payload = json.dumps({"ids": ids}).encode()
    for a in range(retries):
        try:
            req = urllib.request.Request("https://api.gdc.cancer.gov/data", data=payload,
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=600) as r:
                return r.read()
        except Exception as e:
            print(f"  bundle attempt {a+1} failed: {e}", flush=True)
            time.sleep(15 * (a + 1))
    raise RuntimeError("bundle failed after retries")

def parse_maf_bytes(name, b):
    rows = []
    text = b.decode("utf-8", "replace")
    lines = [l for l in text.split("\n") if l and not l.startswith("#")]
    if not lines: return rows
    hdr = lines[0].split("\t")
    col = {c: i for i, c in enumerate(hdr)}
    need = ["Hugo_Symbol", "Variant_Classification", "HGVSp_Short", "Tumor_Sample_Barcode"]
    if not all(c in col for c in need): return rows
    for l in lines[1:]:
        f = l.split("\t")
        if len(f) <= max(col[c] for c in need): continue
        if f[col["Variant_Classification"]] != "Missense_Mutation": continue
        m = HGVSP.match(f[col["HGVSp_Short"]].strip())
        if not m: continue
        rows.append({"sample": f[col["Tumor_Sample_Barcode"]][:16],
                     "gene": f[col["Hugo_Symbol"]],
                     "wt_aa": m.group(1), "pos": int(m.group(2)), "mut_aa": m.group(3),
                     "hgvsp": f[col["HGVSp_Short"]].strip(), "src_file": name})
    return rows

def main():
    man = json.load(open(MAN))
    files = man["files"]
    os.makedirs(OUTDIR, exist_ok=True)
    for bi in range(0, len(files), BUNDLE):
        chunk = files[bi:bi + BUNDLE]
        ck = f"{OUTDIR}/bundle_{bi//BUNDLE:03d}.json"
        if os.path.exists(ck):
            print(f"bundle {bi//BUNDLE} cached", flush=True)
            continue
        ids = [f["file_id"] for f in chunk]
        print(f"bundle {bi//BUNDLE}: {len(ids)} files...", flush=True)
        tar_bytes = post_data(ids)
        muts, seen = [], set()
        with tarfile.open(fileobj=io.BytesIO(tar_bytes)) as tf:
            for member in tf:
                if not member.isfile(): continue
                b = tf.extractfile(member).read()
                if member.name.endswith('.gz'):
                    import gzip as _gz
                    b = _gz.decompress(b)
                muts.extend(parse_maf_bytes(member.name, b))
                seen.add(member.name.split("/")[-1])
        missing = [f["file_name"] for f in chunk if f["file_name"] not in seen]
        json.dump({"bundle": bi // BUNDLE, "file_ids": ids, "n_mutations": len(muts),
                   "missing_files": missing, "mutations": muts}, open(ck, "w"))
        print(f"bundle {bi//BUNDLE}: {len(muts)} missense, {len(missing)} missing", flush=True)
        time.sleep(2)
    print("MAF DOWNLOAD DONE", flush=True)

if __name__ == "__main__":
    main()
