"""Exact intersection of frozen 09d shortlist with HLA Ligand Atlas 2020.12.

Usage: python3 scripts/discovery_hla_atlas_check.py /path/to/hla_2020.12.zip
This comparison is an absence check within a particular benign HLA-ligand release,
not experimental tumor evidence or proof of universal absence.
"""
import csv
import hashlib
import json
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

archive = Path(sys.argv[1]); sha = hashlib.sha256(archive.read_bytes()).hexdigest()
scores = json.load(open('results/discovery/named_candidates_scores.json'))
peptide_scores = defaultdict(list)
for pair, value in scores.items():
    peptide, allele = pair.split('|')
    peptide_scores[peptide].append(float(value))
ranked = sorted(peptide_scores, key=lambda p: (-max(peptide_scores[p]), p))
seen = set(); matches = []
with zipfile.ZipFile(archive) as z, z.open('hla_2020.12/HLA_aggregated.tsv') as raw:
    rows = csv.DictReader((line.decode('utf-8') for line in raw), delimiter='\t')
    for row in rows:
        seen.add(row['peptide_sequence'])
        if row['peptide_sequence'] in peptide_scores:
            matches.append(row['peptide_sequence'])
result = {'source': 'https://hla-ligand-atlas.org/rel/hla_2020.12.zip',
          'archive_sha256': sha, 'release': '2020.12',
          'n_atlas_unique_peptides': len(seen),
          'n_shortlist_unique_peptides': len(ranked),
          'n_top_40_exact_hits': len(set(ranked[:40]) & seen),
          'n_all_exact_hits': len(set(matches)),
          'exact_hits': sorted(set(matches)),
          'interpretation': 'Exact lookup within 2020.12 benign HLA ligand release only; absence is not tumor presentation.'}
out = Path('results/discovery/hla_atlas_2020_12_exact_check.json')
out.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
