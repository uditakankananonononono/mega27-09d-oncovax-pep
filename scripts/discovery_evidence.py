"""Discovery layer: independent evidence for named candidates (addendum
condition 5) + novelty catalog check (condition 1) + GTEx normal-tissue
medians (condition 3 second half). NCBI eutils + PRIDE Archive API +
GTEx API. Politeness-spaced, resume-safe per candidate.
python3 scripts/discovery_evidence.py
"""
import glob, gzip, json, os, time, urllib.parse, urllib.request

D = "results/discovery"
OUT = f"{D}/evidence"

def get(url, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent": "mega27-09d-research/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode()

def pubmed(term):
    q = urllib.parse.quote(term)
    d = json.loads(get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={q}&retmode=json&retmax=5"))
    return d.get("esearchresult", {}).get("idlist", [])

def pride(peptide):
    q = urllib.parse.quote(peptide)
    try:
        d = json.loads(get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword={q}&pageSize=5"))
        hits = d.get("_embedded", {}).get("projects", [])
        return [{"accession": h.get("accession"), "title": (h.get("title") or "")[:160]} for h in hits]
    except Exception as e:
        return {"error": str(e)[:200]}

def gtex(gene):
    try:
        d = json.loads(get(f"https://gtexportal.org/api/v2/reference/gene?geneId={urllib.parse.quote(gene)}"))
        gid = d["data"][0]["gencodeId"] if d.get("data") else None
        if not gid: return {"error": "gene not found"}
        e = json.loads(get(f"https://gtexportal.org/api/v2/expression/medianGeneExpression?gencodeId={urllib.parse.quote(gid)}"))
        tpms = [x.get("median") for x in e.get("data", []) if x.get("median") is not None]
        tpms = sorted(tpms)
        if not tpms: return {"error": "no expression data"}
        return {"gencode_id": gid, "median_across_tissues": tpms[len(tpms)//2],
                "max_across_tissues": tpms[-1], "n_tissues": len(tpms)}
    except Exception as e:
        return {"error": str(e)[:200]}

def main():
    named = json.load(open(f"{D}/named_candidates_scores.json"))
    recs = json.load(gzip.open(f"{D}/funnel_prescr.json.gz", "rt"))["records"]
    # peptide -> carrier records
    carry = {}
    for r in recs:
        for p in r["peptides"]:
            carry.setdefault(p, []).append({"sample": r["sample"], "gene": r["gene"],
                                            "hgvsp": r["hgvsp"], "tpm": r["tpm"]})
    train_peps = set()
    import csv
    for row in csv.DictReader(open("data/processed/splits_v1.csv")):
        train_peps.add(row["peptide"])
    os.makedirs(OUT, exist_ok=True)
    for pair_key, score in sorted(named.items(), key=lambda kv: -kv[1]):
        pep, allele = pair_key.split("|")
        ck = f"{OUT}/{pep}_{allele.replace('*','').replace(':','')}.json"
        if os.path.exists(ck): continue
        ev = {"peptide": pep, "allele": allele, "hgb_score": score,
              "carriers": carry.get(pep, []),
              "novelty": {"in_09d_train_or_test_universe": pep in train_peps}}
        genes = sorted({c["gene"] for c in carry.get(pep, [])})
        terms = [f'"{g}"[Title/Abstract] AND "{c["hgvsp"].split(".")[-1]}"' for c in carry.get(pep, [])[:3] for g in [c["gene"]]]
        ev["pubmed"] = {t: pubmed(t) for t in terms}
        time.sleep(0.4)
        ev["pride"] = pride(pep)
        time.sleep(0.4)
        ev["gtex"] = {g: gtex(g) for g in genes}
        json.dump(ev, open(ck, "w"), indent=1)
        print(f"{pep}|{allele}: pubmed {sum(len(v) for v in ev['pubmed'].values())}, pride {len(ev['pride']) if isinstance(ev['pride'], list) else 'err'}", flush=True)
        time.sleep(0.5)
    print("EVIDENCE DONE", flush=True)

if __name__ == "__main__":
    main()
