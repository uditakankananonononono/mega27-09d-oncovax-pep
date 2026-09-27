#!/usr/bin/env python3
"""Comparator B: DeepImmuno-CNN (frankligy/DeepImmuno, pinned commit
df42ac5b6bddfe531268335e2dcb496559cd488b) on the frozen test 9-10mer
subset. Addendum docs/ADDENDUM_COMPARATOR_B_EXEC_2026-09-28.md.

Phase 1 (infer): reproduce the published architecture + AAindex encoding
verbatim, load the published TF1 checkpoint, score every 9-10mer test
pair. Writes results/comparator_b_deepimmuno_scores.csv +
comparator_b_deepimmuno_meta.json (allele mapping disclosure).
Phase 2 (--score): metrics on the identical 3,890-pair subset for
DeepImmuno, v1 logreg/hgb (from results/model_v1_test_scores.csv) and
NetMHCpan (from results/comparator_a_raw/) with verbatim v1 metric code.
DeepImmuno repo expected at /tmp/deepimmuno (clone path, not committed;
provenance pinned in the addendum and meta JSON).
"""
import csv, json, math, os, sys
import numpy as np
import pandas as pd

DI = "/tmp/deepimmuno"
SPLITS = "data/processed/splits_v1.csv"

# ---- verbatim published architecture (deepimmuno-cnn.py) ----
def seperateCNN():
    import tensorflow as tf
    import tensorflow.keras as keras
    from tensorflow.keras import layers
    input1 = keras.Input(shape=(10, 12, 1))
    input2 = keras.Input(shape=(46, 12, 1))
    x = layers.Conv2D(filters=16, kernel_size=(2, 12))(input1)
    x = layers.BatchNormalization()(x)
    x = keras.activations.relu(x)
    x = layers.Conv2D(filters=32, kernel_size=(2, 1))(x)
    x = layers.BatchNormalization()(x)
    x = keras.activations.relu(x)
    x = layers.MaxPool2D(pool_size=(2, 1), strides=(2, 1))(x)
    x = layers.Flatten()(x)
    x = keras.Model(inputs=input1, outputs=x)
    y = layers.Conv2D(filters=16, kernel_size=(15, 12))(input2)
    y = layers.BatchNormalization()(y)
    y = keras.activations.relu(y)
    y = layers.MaxPool2D(pool_size=(2, 1), strides=(2, 1))(y)
    y = layers.Conv2D(filters=32, kernel_size=(9, 1))(y)
    y = layers.BatchNormalization()(y)
    y = keras.activations.relu(y)
    y = layers.MaxPool2D(pool_size=(2, 1), strides=(2, 1))(y)
    y = layers.Flatten()(y)
    y = keras.Model(inputs=input2, outputs=y)
    combined = layers.concatenate([x.output, y.output])
    z = layers.Dense(128, activation='relu')(combined)
    z = layers.Dropout(0.2)(z)
    z = layers.Dense(1, activation='sigmoid')(z)
    model = keras.Model(inputs=[input1, input2], outputs=z)
    return model

# ---- verbatim published encoding helpers ----
def aaindex(peptide, after_pca):
    amino = 'ARNDCQEGHILKMFPSTWYV-'
    matrix = np.transpose(after_pca)
    encoded = np.empty([len(peptide), 12])
    for i in range(len(peptide)):
        query = peptide[i]
        if query == 'X': query = '-'
        query = query.upper()
        encoded[i, :] = matrix[:, amino.index(query)]
    return encoded

def peptide_data_aaindex(peptide, after_pca):
    length = len(peptide)
    if length == 10:
        encode = aaindex(peptide, after_pca)
    elif length == 9:
        peptide = peptide[:5] + '-' + peptide[5:]
        encode = aaindex(peptide, after_pca)
    else:
        raise ValueError(f"unsupported length {length}")
    encode = encode.reshape(encode.shape[0], encode.shape[1], -1)
    return encode

def dict_inventory(inventory):
    dicA, dicB, dicC = {}, {}, {}
    dic = {'A': dicA, 'B': dicB, 'C': dicC}
    for hla in inventory:
        type_ = hla[4]
        first2 = hla[6:8]
        last2 = hla[8:]
        try:
            dic[type_][first2].append(last2)
        except KeyError:
            dic[type_][first2] = []
            dic[type_][first2].append(last2)
    return dic

def rescue_unknown_hla(hla, dic_inventory):
    type_ = hla[4]
    first2 = hla[6:8]
    last2 = hla[8:]
    big_category = dic_inventory[type_]
    if not big_category.get(first2) == None:
        small_category = big_category.get(first2)
        distance = [abs(int(last2) - int(i)) for i in small_category]
        optimal = min(zip(small_category, distance), key=lambda x: x[1])[0]
        return 'HLA-' + str(type_) + '*' + str(first2) + str(optimal)
    else:
        small_category = list(big_category.keys())
        distance = [abs(int(first2) - int(i)) for i in small_category]
        optimal = min(zip(small_category, distance), key=lambda x: x[1])[0]
        return 'HLA-' + str(type_) + '*' + str(optimal) + str(big_category[optimal][0])

def hla_data_aaindex(hla_dic, hla_type, after_pca, dic_inventory, mapping_log):
    try:
        seq = hla_dic[hla_type]
        mapping_log[hla_type] = hla_type
    except KeyError:
        rescued = rescue_unknown_hla(hla_type, dic_inventory)
        mapping_log[hla_type] = rescued
        seq = hla_dic[rescued]
    encode = aaindex(seq, after_pca)
    encode = encode.reshape(encode.shape[0], encode.shape[1], -1)
    return encode

def infer():
    import tensorflow as tf
    after_pca = np.loadtxt(os.path.join(DI, 'data/after_pca.txt'))
    hla = pd.read_csv(os.path.join(DI, 'data/hla2paratopeTable_aligned.txt'), sep='\t')
    hla_dic = {hla['HLA'].iloc[i]: hla['pseudo'].iloc[i] for i in range(hla.shape[0])}
    inventory = list(hla_dic.keys())
    dic_inv = dict_inventory(inventory)
    model = seperateCNN()
    model.load_weights(os.path.join(DI, 'models/cnn_model_331_3_7/'))

    rows = list(csv.DictReader(open(SPLITS)))
    test = [r for r in rows if r['split'] == 'test' and len(r['peptide']) in (9, 10)]
    pairs = sorted({(r['peptide'], r['allele']) for r in test})
    print(f"infer: {len(pairs)} unique 9-10mer pairs", flush=True)

    mapping_log = {}
    out_rows = []
    B = 512
    for i in range(0, len(pairs), B):
        chunk = pairs[i:i + B]
        x1 = np.stack([peptide_data_aaindex(p, after_pca) for p, a in chunk])
        x2 = np.stack([hla_data_aaindex(hla_dic, a.replace(':', ''), after_pca, dic_inv, mapping_log) for p, a in chunk])
        scores = model.predict(x=[x1, x2], verbose=0).reshape(-1)
        for (p, a), s in zip(chunk, scores):
            out_rows.append((p, a, float(s)))
        print(f"infer {min(i + B, len(pairs))}/{len(pairs)}", flush=True)

    with open('results/comparator_b_deepimmuno_scores.csv', 'w') as fh:
        fh.write('peptide,allele,deepimmuno_score\n')
        for p, a, s in out_rows:
            fh.write(f"{p},{a},{s}\n")
    norm_map = {}
    for k, v in mapping_log.items():
        pretty = k[:6] + ':' + k[6:] if len(k) > 6 else k
        pretty_v = v[:6] + ':' + v[6:] if len(v) > 6 else v
        norm_map[pretty] = pretty_v
    json.dump({"allele_mapping": norm_map,
               "deepimmuno_commit": "df42ac5b6bddfe531268335e2dcb496559cd488b",
               "weights": "models/cnn_model_331_3_7 (TF1 checkpoint)",
               "n_pairs": len(out_rows)},
              open('results/comparator_b_deepimmuno_meta.json', 'w'), indent=1)
    print("INFER DONE", flush=True)

def ndcg_at_k(score, gain, k=100):
    order = np.argsort(-score)[:k]
    dcg = sum((2**gain[i] - 1) / math.log2(r + 2) for r, i in enumerate(order))
    ideal = np.argsort(-gain)[:k]
    idcg = sum((2**gain[i] - 1) / math.log2(r + 2) for r, i in enumerate(ideal))
    return dcg / idcg if idcg > 0 else 0.0

def metric_block(s, y, g):
    from sklearn.metrics import average_precision_score, roc_auc_score
    k = min(100, len(s))
    order = np.argsort(-s)[:k]
    return dict(AUPRC=float(average_precision_score(y, s)),
                AUROC=float(roc_auc_score(y, s)),
                Recall_at_100=float(y[order].mean()),
                NDCG_at_100=float(ndcg_at_k(s, g, k)))

def score():
    rows = list(csv.DictReader(open(SPLITS)))
    test = [r for r in rows if r['split'] == 'test' and len(r['peptide']) in (9, 10)]
    di = {(r['peptide'], r['allele']): float(r['deepimmuno_score'])
          for r in csv.DictReader(open('results/comparator_b_deepimmuno_scores.csv'))}
    v1 = {}
    for r in csv.DictReader(open('results/model_v1_test_scores.csv')):
        v1[(r['peptide'], r['allele'])] = (float(r['logreg_score']), float(r['hgb_score']))
    nm = {}
    for fn in os.listdir('results/comparator_a_raw'):
        if not fn.endswith('.json'):
            continue
        d = json.load(open(f'results/comparator_a_raw/{fn}'))
        if d.get('unsupported'):
            continue
        for t in d['data']['results']:
            if t['type'] != 'peptide_table':
                continue
            cols = [c['name'] for c in t['table_columns']]
            ip, ia = cols.index('peptide'), cols.index('allele')
            ipct = cols.index('netmhcpan_el_percentile')
            for row in t['table_data']:
                nm[(row[ip], row[ia])] = float(row[ipct])
    kept, missing = [], []
    for r in test:
        key = (r['peptide'], r['allele'])
        if key in di and key in v1 and key in nm:
            kept.append(r)
        else:
            missing.append(key)
    print(f"subset pairs scored by all systems: {len(kept)}; missing: {len(missing)}", flush=True)
    if missing:
        print("sample missing:", missing[:5], flush=True)
    y = np.array([int(r['label']) for r in kept])
    g = np.array([float(r['positive_fraction']) for r in kept])
    keys = [(r['peptide'], r['allele']) for r in kept]
    out = {
        "deepimmuno": metric_block(np.array([di[k] for k in keys]), y, g),
        "netmhcpan_el_rank": metric_block(np.array([-nm[k] for k in keys]), y, g),
        "v1_logreg": metric_block(np.array([v1[k][0] for k in keys]), y, g),
        "v1_hgb": metric_block(np.array([v1[k][1] for k in keys]), y, g),
    }
    res = dict(comparator="DeepImmuno-CNN (pinned commit df42ac5b)",
               addendum="docs/ADDENDUM_COMPARATOR_B_EXEC_2026-09-28.md",
               n_test=len(kept), n_missing_any_system=len(missing),
               test_pos_rate=float(y.mean()), results=out)
    json.dump(res, open('results/comparator_b_deepimmuno.json', 'w'), indent=1)
    print(json.dumps(out, indent=1), flush=True)
    print("SCORE DONE", flush=True)

if __name__ == "__main__":
    if "--score" in sys.argv:
        score()
    else:
        infer()
