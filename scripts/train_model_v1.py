"""09d model v1: immunogenicity predictor trained on frozen train split only.
Features: MHCflurry-2.2.1 predicted affinity (log nM) for the exact allele,
AA composition (20), length, Kyte-Doolittle hydrophobicity mean/range, net charge,
aromatic fraction. Models: logistic regression (baseline) + HistGradientBoosting.
Locked metrics on frozen test: AUPRC (primary), Recall@100, NDCG@100 (graded rel =
positive_fraction). NO test-set peeking during training; single evaluation at end.
"""
import csv, json, math, pickle
import numpy as np
from collections import Counter

KD = dict(A=1.8,C=2.5,D=-3.5,E=-3.5,F=2.8,G=-0.4,H=-3.2,I=4.5,K=-3.9,L=3.8,
          M=1.9,N=-3.5,P=-1.6,Q=-3.5,R=-4.5,S=-0.8,T=-0.7,V=4.2,W=-0.9,Y=-1.3)
CH = dict(D=-1,E=-1,K=1,R=1,H=0.1)
ARO = set('FWY')
AAS = 'ACDEFGHIKLMNPQRSTVWY'

rows = list(csv.DictReader(open('data/processed/splits_v1.csv')))
tr = [r for r in rows if r['split'] == 'train']
te = [r for r in rows if r['split'] == 'test']

# --- MHCflurry affinity for all needed pairs (cached) ---
from mhcflurry import Class1AffinityPredictor
pred = Class1AffinityPredictor.load()
pairs = sorted({(r['peptide'], r['allele']) for r in tr + te})
peps = sorted({p for p, a in pairs})
alleles = sorted({a for p, a in pairs})
need = {}
res = pred.predict_to_dataframe(peptides=[p for p, a in pairs],
                                alleles=[a for p, a in pairs], throw=False)
for (p, a), aff in zip([(r_.peptide, r_.allele) for r_ in res.itertuples()], res.prediction):
    need[(p, a)] = aff
pickle.dump(need, open('/tmp/09d_aff_cache.pkl', 'wb'))
print('affinities done', len(need), flush=True)

def feats(r):
    s = r['peptide']; n = len(s)
    c = Counter(s)
    comp = [c.get(a, 0) / n for a in AAS]
    kd = [KD[x] for x in s]
    aff = need.get((s, r['allele']), 500.0)
    return comp + [n, float(np.mean(kd)), float(np.max(kd) - np.min(kd)),
                   sum(CH.get(x, 0) for x in s),
                   sum(1 for x in s if x in ARO) / n,
                   math.log(max(aff, 1e-3))]

Xtr = np.array([feats(r) for r in tr]); ytr = np.array([int(r['label']) for r in tr])
Xte = np.array([feats(r) for r in te]); yte = np.array([int(r['label']) for r in te])
gte = np.array([float(r['positive_fraction']) for r in te])
pickle.dump(dict(Xtr=Xtr, ytr=ytr, Xte=Xte, yte=yte, gte=gte), open('/tmp/09d_v1_ckpt.pkl', 'wb'))
print('features done', Xtr.shape, Xte.shape, flush=True)

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score

def ndcg_at_k(score, gain, k=100):
    order = np.argsort(-score)[:k]
    dcg = sum((2**gain[i] - 1) / math.log2(r + 2) for r, i in enumerate(order))
    ideal = np.argsort(-gain)[:k]
    idcg = sum((2**gain[i] - 1) / math.log2(r + 2) for r, i in enumerate(ideal))
    return dcg / idcg if idcg > 0 else 0.0

out = {}
for name, clf in [('logreg', LogisticRegression(max_iter=2000, C=1.0)),
                  ('hgb', HistGradientBoostingClassifier(max_iter=300, learning_rate=0.08,
                                                         max_depth=None, random_state=2709))]:
    clf.fit(Xtr, ytr)
    s = clf.predict_proba(Xte)[:, 1]
    k = min(100, len(s))
    order = np.argsort(-s)[:k]
    out[name] = dict(AUPRC=float(average_precision_score(yte, s)),
                     AUROC=float(roc_auc_score(yte, s)),
                     Recall_at_100=float(yte[order].mean()),
                     NDCG_at_100=float(ndcg_at_k(s, gte, k)))
    print(name, out[name], flush=True)
    pickle.dump(clf, open(f'/tmp/09d_model_{name}.pkl', 'wb'))

json.dump(dict(n_train=len(tr), n_test=len(te),
               test_pos_rate=float(yte.mean()), results=out),
          open('results/model_v1_frozen_test.json', 'w'), indent=1)
print('DONE', flush=True)
