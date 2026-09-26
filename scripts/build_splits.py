"""Locked split builder for 09d: allele-group + peptide-cluster holdout.
Prereg: no peptide cluster or allele family shared between train/test. Frozen before training.
Peptide clusters: single-linkage over shared exact 7-mers (union-find), k=7 chosen as
specific enough to avoid transitive lumping while catching near-duplicates.
Allele families: two-field group root e.g. A*02, B*07 (HLA-A*02:01 -> A*02).
Partition rule: TEST iff cluster_id in test_clusters AND family in test_families;
TRAIN iff cluster in train_clusters AND family in train_families; else UNASSIGNED (excluded).
Seed fixed. This script runs BEFORE any model training; output committed as frozen artifact.
"""
import csv, json, random, re
from collections import defaultdict

random.seed(2709)
K = 7
TEST_FRAC_CLUSTERS = 0.20
TEST_FRAC_FAMILIES = 0.20

rows = list(csv.DictReader(open('data/processed/iedb_tcell_class1_human.csv')))
peps = sorted({r['peptide'] for r in rows})

# union-find over peptides via shared k-mers
parent = {p: p for p in peps}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[ra] = rb

idx = defaultdict(list)
for p in peps:
    for i in range(len(p) - K + 1):
        idx[p[i:i+K]].append(p)
for ps in idx.values():
    for q in ps[1:]:
        union(ps[0], q)

clusters = defaultdict(list)
for p in peps:
    clusters[find(p)].append(p)
cids = sorted(clusters, key=lambda c: (-len(clusters[c]), c))
cluster_of = {p: f'C{rank:05d}' for rank, c in enumerate(cids) for p in clusters[c]}

def fam(a):
    m = re.match(r'^HLA-([A-C])\*(\d{2})', a)
    return f'{m.group(1)}*{m.group(2)}'
fams = sorted({fam(r['allele']) for r in rows})

# positive-aware cluster split: stratify by cluster positive mass
_cs = defaultdict(int); _cn = defaultdict(int)
for r in rows:
    c = cluster_of[r['peptide']]
    _cs[c] += int(r['label']); _cn[c] += 1
cmass = [(cluster_of[clusters[c][0]], _cs[cluster_of[clusters[c][0]]], _cn[cluster_of[clusters[c][0]]]) for c in cids]
# stratified interleaved deal by positive rate: deterministic, keeps test pos-rate near train
cmass.sort(key=lambda t: (-(t[1] + 0.5) / (t[2] + 1), t[0]))
test_clusters, train_clusters = set(), set()
for i, (c, s, n) in enumerate(cmass):
    (test_clusters if i % 5 == 0 else train_clusters).add(c)

random.shuffle(fams)
ntf = max(1, round(TEST_FRAC_FAMILIES * len(fams)))
test_families, train_families = set(fams[:ntf]), set(fams[ntf:])

n_test = n_train = n_un = 0
with open('data/processed/splits_v1.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['peptide', 'allele', 'label', 'positive_fraction', 'n_assays',
                'cluster_id', 'allele_family', 'split'])
    for r in rows:
        cid = cluster_of[r['peptide']]; fa = fam(r['allele'])
        if cid in test_clusters and fa in test_families:
            sp = 'test'; n_test += 1
        elif cid in train_clusters and fa in train_families:
            sp = 'train'; n_train += 1
        else:
            sp = 'unassigned'; n_un += 1
        w.writerow([r['peptide'], r['allele'], r['label'], r['positive_fraction'],
                    r['n_assays'], cid, fa, sp])

spec = dict(seed=2709, kmer=K, n_peptides=len(peps), n_clusters=len(cids),
            n_test_clusters=len(test_clusters), n_families=len(fams),
            test_families=sorted(test_families), train_families=sorted(train_families),
            n_train=n_train, n_test=n_test, n_unassigned=n_un,
            rule='TEST iff test cluster AND test family; TRAIN iff train cluster AND train family; else unassigned/excluded')
json.dump(spec, open('data/processed/splits_v1_spec.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in spec.items() if not isinstance(v, list)}, indent=1))
