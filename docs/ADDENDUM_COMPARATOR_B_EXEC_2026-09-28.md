# Addendum: Comparator B execution mechanics (2026-09-28, locked before scoring)

The 2026-09-26 prereg locked Comparator B (a published immunogenicity
predictor - DeepImmuno, PMC7781330 - on the identical test partition) and
the metrics (AUPRC primary, Recall@100, NDCG@100). This addendum locks
ONLY the execution mechanics, written and committed BEFORE any DeepImmuno
score is compared to labels. Feasibility checks so far touched code and
tables only - no test labels were scored.

## Tool (locked)
- DeepImmuno-CNN as published: github.com/frankligy/DeepImmuno, commit
  pinned at clone time (recorded in the runner output), model weights
  models/cnn_model_331_3_7 (TF1 checkpoint), AAindex1/after_pca encoding
  and the published hla2paratopeTable_aligned.txt pseudo-sequence table.
  Inference uses the published architecture and encoding verbatim.
- Runtime: tensorflow-cpu via pip on the sandbox (published env is
  TF 2.3 / py3.6); any compatibility shims needed to load the checkpoint
  under the installed TF version will be enumerated in the results JSON,
  never silently patched.

## Scope discovered during feasibility (disclosed, locked)
- The published peptide encoder (peptide_data_aaindex) supports EXACTLY
  9- and 10-mer peptides; other lengths are outside the published model's
  input class. DeepImmuno is therefore scored on the 3,890 9-10-mer test
  pairs. The 459 pairs of other lengths (70 x 8mer, 143 x 11mer, 14 x
  12mer, 8 x 13mer, 189 x 14mer, 35 x 15mer) are excluded SYMMETRICALLY:
  the Comparator-B head-to-head recomputes model v1 (logreg, hgb) AND
  NetMHCpan metrics on the identical 3,890-pair subset from the committed
  per-row scores (results/model_v1_test_scores.csv and
  results/comparator_a_raw/), using the verbatim v1 metric code.
- HLA handling: colon-normalized exact match against the published
  pseudo-sequence table (13/32 test alleles); the remaining 19 alleles
  are mapped by the published rescue_unknown_hla function (nearest
  allele in the table). The full allele mapping is disclosed in the
  results JSON.
- Comparator partitions honestly differ by tool: Comparator A (NetMHCpan)
  covers 8-14-mers (n=4,314), Comparator B (DeepImmuno) covers 9-10-mers
  (n=3,890). Each head-to-head uses that comparator's identical subset.
  The prereg BEAT rule is applied per comparator on its own partition
  and reported as such - no cross-partition claims.

## Decision rule (re-stated from prereg, unchanged)
BEAT requires exceeding BOTH comparators on AUPRC on the identical frozen
partition for that comparator. Outcomes reported as-is either way;
negatives are not terminal. If DeepImmuno cannot be executed faithfully,
NeoTImmuML (PMC12585993) is the documented fallback under a new addendum.

## Supplement 1 (2026-09-28, pre-execution, label-free): Keras-3 weight-load shim

The published weights (`models/cnn_model_331_3_7`) are a TF2 object-graph
checkpoint (variables under generic `layer_with_weights-N` paths; files
`.data-00000-of-00001`/`.index`). tensorflow-cpu 2.21.0 (Keras 3) cannot
load it natively. The following shim is locked here BEFORE any execution
and before any label contact:

1. File-level repack only: copy the two published checkpoint files to
   /tmp/di_ckpt/ckpt.* with a rewritten `checkpoint` pointer file.
   Bytes unchanged; no numeric modification of any weight.
2. Conv and Dense layers are assigned unambiguously by kernel shape:
   [15,12,1,16]=HLA conv1 (lww-0), [2,12,1,16]=peptide conv1 (lww-1),
   [2,1,16,32]=peptide conv2 (lww-4), [9,1,16,32]=HLA conv2 (lww-5),
   [256,128]=dense (lww-8), [128,1]=dense_1 (lww-9). Shape asserts on
   every assignment.
3. The four BatchNormalization layers are shape-ambiguous in pairs
   (16,16) lww-2/3 and (32,32) lww-6/7. Pairing is resolved label-free:
   empirical pre-BN channel means computed on unlabeled probe encodings
   (first 1024 pairs' input features only; labels are never read in the
   shim) are matched to each candidate BN's stored
   moving_mean/moving_variance by z-scored MSE. BN16 pair is resolved
   first on first-conv outputs (BN-independent); BN32 pair is then
   resolved on second-conv outputs with BN16 fixed. Assignment requires
   a decisive margin (cost ratio >= 2) or the run aborts. Assignment,
   distances and margins are disclosed in
   results/comparator_b_deepimmuno_meta.json.
4. Sanity gate: probe predictions must spread within (0,1) with std >
   1e-4; spread recorded in the meta JSON. Saturation aborts the run.

Comparator metrics are computed only afterwards by the locked --score
phase from committed per-row scores.
