"""Exact shortlist overlap with dbPepNeo 2019 HC/MC downloadable catalogs.
Usage: python3 scripts/discovery_dbpepneo_check.py /path/to/HC.zip /path/to/MC.zip
No-hit is bounded to these small historical catalogs; it is not discovery evidence.
"""
import hashlib, io, json, sys, zipfile
from pathlib import Path
import openpyxl

shortlist = {x.split('|')[0] for x in json.load(open('results/discovery/named_candidates_scores.json'))}
result = {'source': 'http://www.biostatistics.online/dbPepNeo/download.html',
          'n_shortlist_unique_peptides': len(shortlist), 'catalogs': {}}
for label, path in zip(('HC', 'MC'), sys.argv[1:]):
    p = Path(path)
    with zipfile.ZipFile(p) as z:
        names = [n for n in z.namelist() if n.endswith('.xlsx')]
        assert len(names) == 1
        wb = openpyxl.load_workbook(io.BytesIO(z.read(names[0])), read_only=True, data_only=True)
        sheet = wb.active
        rows = list(sheet.values)[2:]
    peptides = {str(r[3]).strip().upper() for r in rows if r[3]}
    result['catalogs'][label] = {'archive_sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
                                  'n_data_rows': len(rows), 'n_unique_peptides': len(peptides),
                                  'n_shortlist_exact_hits': len(peptides & shortlist),
                                  'exact_hits': sorted(peptides & shortlist)}
result['interpretation'] = 'Exact lookup in dbPepNeo 2019 catalogs only; absence does not prove novelty or HLA presentation.'
Path('results/discovery/dbpepneo_2019_exact_check.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
