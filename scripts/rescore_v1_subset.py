#!/usr/bin/env python3
"""09d model v1 RESCORE on the symmetric-exclusion subset (8-14-mers, n=4314).
Mirrors scripts/train_model_v1.py exactly (same features, models, seed,
metric code - copied verbatim, disclosed) with three additions:
  1. per-row test scores dumped to results/model_v1_test_scores.csv
  2. metrics on the full test set (reproducibility check vs committed
     model_v1_frozen_test.json) AND on the 8-14-mer subset (head-to-head
     basis for Comparator A, addendum note 1)
  3. affinity cache kept at /tmp/09d_aff_cache.pkl (reuse if present)
Run AFTER mhcflurry is installed. Single evaluation; no test peeking.
"""
import csv, json, math, os, pickle
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

from mhcflurry import Class1AffinityPredictor
if os.path.exists('/tmp/09d_aff_cache.pkl'):
    need = pickle.load(open('/tmp/09d_aff_cache.pkl', 'rb'))
    print('affinity cache loaded', len(need), flush=True)
else:
    pred = Class1AffinityPredictor.load()
    pairs = sorted({(r['peptide'], r['allele']) for r in tr + te})
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

def metric_block(s, y, g):
    k = min(100, len(s))
    order = np.argsort(-s)[:k]
    return dict(AUPRC=float(average_precision_score(y, s)),
                AUROC=float(roc_auc_score(y, s)),
                Recall_at_100=float(y[order].mean()),
                NDCG_at_100=float(ndcg_at_k(s, g, k)))

sub = np.array([len(r['peptide']) <= 14 for r in te])  # symmetric-exclusion subset
out = {}
with open('results/model_v1_test_scores.csv', 'w') as fh:
    fh.write('peptide,allele,label,positive_fraction,logreg_score,hgb_score\n')
    score_cols = {}
    for name, clf in [('logreg', LogisticRegression(max_iter=2000, C=1.0)),
                      ('hgb', HistGradientBoostingClassifier(max_iter=300, learning_rate=0.08,
                                                             max_depth=None, random_state=2709))]:
        clf.fit(Xtr, ytr)
        s = clf.predict_proba(Xte)[:, 1]
        score_cols[name] = s
        out[name] = dict(full_test=metric_block(s, yte, gte),
                         subset_8_14=metric_block(s[sub], yte[sub], gte[sub]))
        print(name, json.dumps(out[name]), flush=True)
    for i, r in enumerate(te):
        fh.write(f"{r['peptide']},{r['allele']},{r['label']},{r['positive_fraction']},"
                 f"{score_cols['logreg'][i]},{score_cols['hgb'][i]}\n")

json.dump(dict(note="rescore for Comparator A head-to-head subset; mirrors train_model_v1.py",
               n_train=len(tr), n_test=len(te), n_subset=int(sub.sum()),
               results=out),
          open('results/model_v1_rescore_subset.json', 'w'), indent=1)
print('RESCORE DONE', flush=True)
