import csv, json
from collections import defaultdict
AA = set("ACDEFGHIKLMNPQRSTVWY")
rows = defaultdict(list)
n_total = 0
with open('data/raw/tcell_full_v3.csv') as f:
    r = csv.reader(f)
    next(r); next(r)
    for row in r:
        n_total += 1
        if n_total % 500000 == 0:
            print(f"scanned {n_total}, kept {sum(len(v) for v in rows.values())}", flush=True)
        try:
            pep, host, qual, allele, mclass = row[11], row[43], row[122], row[141], row[145]
        except IndexError:
            continue
        if host != 'Homo sapiens (human)' or mclass != 'I' or not __import__('re').match(r'^HLA-[A-C]\*\d{2}:\d{2}[A-Z]?$', allele):
            continue
        if not (8 <= len(pep) <= 15) or not set(pep) <= AA:
            continue
        if qual.startswith('Positive'):
            lab = 1
        elif qual == 'Negative':
            lab = 0
        else:
            continue
        rows[(pep, allele)].append(lab)
print(f"total scanned {n_total}", flush=True)
out = [(p, a, round(sum(v)/len(v), 3), len(v), 1 if sum(v)/len(v) >= 0.5 else 0)
       for (p, a), v in rows.items()]
with open('data/processed/iedb_tcell_class1_human.csv', 'w') as f:
    f.write('peptide,allele,positive_fraction,n_assays,label\n')
    for p, a, frac, n, lab in sorted(out):
        f.write(f'{p},{a},{frac},{n},{lab}\n')
stats = {"pairs": len(out), "positives": sum(1 for *_, l in out if l == 1),
         "alleles": len({a for _, a, *_ in out})}
json.dump(stats, open('data/processed/build_stats.json', 'w'), indent=1)
print(stats, flush=True)
