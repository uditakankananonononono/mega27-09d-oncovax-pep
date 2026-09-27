# Addendum: Comparator A execution mechanics (2026-09-28, locked before scoring)

The 2026-09-26 prereg locked Comparator A (NetMHCpan-4.1 binding rank on
the identical frozen test set) and the metrics (AUPRC, Recall@100,
NDCG@100). This addendum locks ONLY the execution mechanics, before any
comparator score is compared to labels. The earlier IEDB-API 403 egress
block cleared 2026-09-28 (live-verified); the prereg's "run via free
IEDB tools" route is now executable, so no proxy substitution is needed.

## Mechanics (locked)
- Source: IEDB nextgen tools API, POST /api/v1/pipeline, tool_group mhci,
  predictor type=binding method=netmhcpan_el (NetMHCpan-4.1 EL, the
  current IEDB default for this method name; disclosed as EL, not BA).
- Batching: per test allele (32) x exact peptide length L in {8..15};
  each peptide submitted as its own FASTA record with
  peptide_length_range [L,L], so every record yields exactly itself and
  no spurious k-mer windows enter any score.
- Score per (peptide, allele) test pair: the API's percentile rank for
  the exact peptide string (lower rank = stronger binder). Model-side
  comparator direction: immunogenicity score = -rank (so larger =
  better), and AUPRC/Recall@100/NDCG@100 computed with the same code
  path as model v1 (scripts/train_model_v1.py metric functions or a
  verbatim copy disclosed in the runner).
- Missing API allele: if any of the 32 test alleles is rejected by the
  API, the affected pairs are excluded from BOTH comparator and model
  metrics for the head-to-head (disclosed count), never silently.
- Retries: transient API failures retried up to 3x with backoff; pairs
  still failing after 3 rounds are disclosed and excluded symmetrically.
- Raw API responses cached under results/comparator_a_raw/ (committed
  unless size-prohibitive; then hashes + regeneration script committed).

## Decision rule (re-stated from prereg, unchanged)
BEAT requires exceeding BOTH comparators on AUPRC on the identical
frozen test set. Comparator A alone passing/failing is reported as-is;
Comparator B (DeepImmuno) follows under its own execution addendum.

## Discovery note 1 (2026-09-28, before any scoring): 15-mers unsupported
Live API response: "netmhcpan cannot predict binding for allele HLA-A*01:01
with all requested peptide lengths" - NetMHCpan-4.1 EL via IEDB supports
8-14-mers only. The 35 L15 test pairs (27 positive, disclosed: 77% positive
rate, so exclusion is not label-neutral) are excluded SYMMETRICALLY from
comparator and model metrics per this addendum's exclusion rule. The
head-to-head set is therefore 4,314 pairs (all 8-14-mers), and model v1
metrics will be recomputed on exactly that subset (per-row rescore,
same seed/features) rather than quoted from the full-set run. The full-set
v1 numbers remain committed as the full-set reference.
