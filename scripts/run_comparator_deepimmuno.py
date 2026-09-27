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

# ---- Keras-3 compat shim (addendum supplement 1; label-free) ----
CKPT_REPACK = "/tmp/di_ckpt"

def _repack_checkpoint():
    import shutil
    srcd = os.path.join(DI, 'models/cnn_model_331_3_7')
    os.makedirs(CKPT_REPACK, exist_ok=True)
    for a, b in [('.data-00000-of-00001', 'ckpt.data-00000-of-00001'),
                 ('.index', 'ckpt.index')]:
        shutil.copyfile(os.path.join(srcd, a), os.path.join(CKPT_REPACK, b))
    with open(os.path.join(CKPT_REPACK, 'checkpoint'), 'w') as fh:
        fh.write('model_checkpoint_path: "ckpt"\nall_model_checkpoint_paths: "ckpt"\n')

def load_published_weights(model, x1_probe, x2_probe):
    import tensorflow as tf
    _repack_checkpoint()
    reader = tf.train.load_checkpoint(os.path.join(CKPT_REPACK, 'ckpt'))
    def var(i, name):
        return reader.get_tensor(f'layer_with_weights-{i}/{name}/.ATTRIBUTES/VARIABLE_VALUE')
    def assign(layer_name, i):
        layer = model.get_layer(layer_name)
        got = [var(i, 'kernel'), var(i, 'bias')]
        for w, g in zip(layer.get_weights(), got):
            assert w.shape == g.shape, f"{layer_name}: checkpoint shape {g.shape} != model {w.shape}"
        layer.set_weights(got)
    # convs + head unambiguous by kernel shape (addendum supplement 1)
    assign('conv2d', 1)      # peptide conv1 [2,12,1,16]
    assign('conv2d_2', 0)    # HLA conv1 [15,12,1,16]
    assign('conv2d_1', 4)    # peptide conv2 [2,1,16,32]
    assign('conv2d_3', 5)    # HLA conv2 [9,1,16,32]
    assign('dense', 8)       # [256,128]
    assign('dense_1', 9)     # [128,1]
    report = {"shim": "Keras3 cannot read TF2 object checkpoint; file-level repack to /tmp/di_ckpt (bytes unchanged, prefix renamed) + name-mapped variable assignment; no numeric modification of weights",
              "conv_dense_assignment": {"conv2d(pep1)": 1, "conv2d_2(hla1)": 0, "conv2d_1(pep2)": 4, "conv2d_3(hla2)": 5, "dense": 8, "dense_1": 9},
              "bn_pairing_method": "label-free: empirical pre-BN channel means on unlabeled probe encodings matched to stored moving_mean/moving_variance via z-MSE; aborts if margin < 2"}
    def chmean(t):
        return np.asarray(t).mean(axis=(0, 1, 2))
    def zdist(m, i):
        mu = var(i, 'moving_mean'); v = var(i, 'moving_variance')
        return float(((m - mu) ** 2 / (v + 1e-5)).mean())
    def pair(pep_layer, hla_layer, m_pep, m_hla, cand):
        d = {i: (zdist(m_pep, i), zdist(m_hla, i)) for i in cand}
        costA = d[cand[0]][0] + d[cand[1]][1]
        costB = d[cand[1]][0] + d[cand[0]][1]
        pep_i, hla_i = (cand[0], cand[1]) if costA <= costB else (cand[1], cand[0])
        margin = max(costA, costB) / (min(costA, costB) + 1e-12)
        assert margin >= 2.0, f"BN pairing indecisive: margin {margin}"
        for layer_name, i in [(pep_layer, pep_i), (hla_layer, hla_i)]:
            model.get_layer(layer_name).set_weights(
                [var(i, 'gamma'), var(i, 'beta'), var(i, 'moving_mean'), var(i, 'moving_variance')])
        return {"pep": pep_i, "hla": hla_i, "margin": margin,
                "zdist": {str(k): [round(x, 6) for x in v] for k, v in d.items()}}
    L = model.get_layer
    # BN16 pair: first-conv outputs are BN-independent
    m_pep1 = chmean(L('conv2d')(x1_probe))
    m_hla1 = chmean(L('conv2d_2')(x2_probe))
    report['bn16_assignment'] = pair('batch_normalization', 'batch_normalization_2', m_pep1, m_hla1, (2, 3))
    # BN32 pair: second-conv outputs with BN16 fixed
    bn = lambda name, t: L(name)(t, training=False)
    o_pep = L('conv2d_1')(tf.nn.relu(bn('batch_normalization', L('conv2d')(x1_probe))))
    pool = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=(2, 1))
    o_hla = L('conv2d_3')(pool(tf.nn.relu(bn('batch_normalization_2', L('conv2d_2')(x2_probe)))))
    report['bn32_assignment'] = pair('batch_normalization_1', 'batch_normalization_3', chmean(o_pep), chmean(o_hla), (6, 7))
    report['tf_version'] = tf.__version__
    return report

def infer():
    import tensorflow as tf
    after_pca = np.loadtxt(os.path.join(DI, 'data/after_pca.txt'))
    hla = pd.read_csv(os.path.join(DI, 'data/hla2paratopeTable_aligned.txt'), sep='\t')
    hla_dic = {hla['HLA'].iloc[i]: hla['pseudo'].iloc[i] for i in range(hla.shape[0])}
    inventory = list(hla_dic.keys())
    dic_inv = dict_inventory(inventory)
    rows = list(csv.DictReader(open(SPLITS)))
    test = [r for r in rows if r['split'] == 'test' and len(r['peptide']) in (9, 10)]
    pairs = sorted({(r['peptide'], r['allele']) for r in test})
    print(f"infer: {len(pairs)} unique 9-10mer pairs", flush=True)

    model = seperateCNN()
    probe = pairs[:1024]
    x1p = np.stack([peptide_data_aaindex(p, after_pca) for p, a in probe])
    x2p = np.stack([hla_data_aaindex(hla_dic, a.replace(':', ''), after_pca, dic_inv, {}) for p, a in probe])
    shim_report = load_published_weights(model, x1p, x2p)
    probe_pred = model.predict(x=[x1p, x2p], verbose=0).reshape(-1)
    shim_report['probe_score_spread'] = {"min": float(probe_pred.min()), "max": float(probe_pred.max()),
                                         "std": float(probe_pred.std())}
    assert probe_pred.std() > 1e-4 and 0 <= probe_pred.min() and probe_pred.max() <= 1, "saturated or degenerate predictions after shim"
    print("shim:", json.dumps(shim_report), flush=True)

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
    def pretty(k):
        return k[:6] + k[6:8] + ':' + k[8:] if len(k) > 8 else k
    norm_map = {pretty(k): pretty(v) for k, v in mapping_log.items()}
    json.dump({"allele_mapping": norm_map,
               "deepimmuno_commit": "df42ac5b6bddfe531268335e2dcb496559cd488b",
               "weights": "models/cnn_model_331_3_7 (TF2 object checkpoint df42ac5)",
               "keras3_shim": shim_report,
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
