# ADDENDUM - DISCOVERY LAYER EXECUTION (2026-09-28, locked before any candidate scoring)

Parent direction: main-agent message 2026-09-28 03:52 IST ("start the 09d
discovery layer (TCGA/candidate generation) after C6 combines. Approved.").
Executes the locked discovery criterion in docs/PREREGISTRATION.md
(5 conditions, unchanged). Feasibility pass (label-free) completed before
this addendum: GDC open MAF + STAR-count rails verified reachable from the
sandbox; Xena hubs 403 from sandbox (disclosed); GTEx API reachable.

## Cohort (locked)
TCGA-SKCM (cutaneous melanoma): highest somatic mutation burden among
common TCGA cohorts (public knowledge), 472 open masked-MAF aliquots and
473 open STAR gene-count files (GDC files API, counted 2026-09-28).
Single-cohort design; no cohort cherry-picking after outcomes.

## Data sources (locked, all open-access)
- Somatic mutations: GDC per-aliquot Masked Somatic Mutation MAFs
  (annotation spec gdc-2.0.0-aliquot-merged-masked), project TCGA-SKCM,
  via api.gdc.cancer.gov. File manifest (UUID list + md5) committed.
- Tumor expression: GDC per-aliquot STAR gene counts (same project).
  Streamed per file; only mutated-gene rows retained; raw files deleted
  after extraction (box has 2GB RAM / limited disk). Provenance: file
  UUIDs + md5 recorded per sample.
- Normal-tissue expression: GTEx portal API v2 (gtexportal.org),
  per-gene median TPM across normal tissues.
- Self-proteome filter: UniProt reviewed human proteome (exact source +
  download date recorded at execution; if unreachable, disclosed and the
  filter falls back to the IEDB-hosted reference or the run pauses).
- Novelty catalogs (condition 1): exact-match against (a) the IEDB-derived
  peptide set already committed in this repo (data/processed/splits_v1.csv
  universe + training peptides) and (b) the full IEDB T-cell assay export
  if reachable at execution; each catalog checked is disclosed with its
  lock timestamp. No post-hoc catalog additions.

## Funnel (locked thresholds, applied in order)
1. Mutations: Variant_Classification == Missense_Mutation with a valid
   HGVSp single-residue substitution. Indels/frameshifts/silent EXCLUDED
   this round (disclosed; model is substitution-window based).
2. Neopeptides: all 8-15mer substrings of the +/-7 window around the
   mutated residue that contain the mutated residue; dedup by
   (peptide, allele-panel) downstream; wild-type-mutant pairs recorded.
3. Expression pre-filter: mutated gene >= 1 TPM in the mutation-carrying
   tumor sample (STAR counts, locked threshold).
4. Self filter: peptide must not exactly match any 8-15mer substring of
   the reviewed human reference proteome (removes non-mutated sequences).
5. Presentation pre-screen: MHCflurry-2.2.1 affinity percentile <= 2.0
   for at least one panel allele (panel locked below). Engine + version
   disclosed; this is a pre-screen, not the comparator.
6. Model score: v1 hgb (primary) + v1 logreg (secondary, reported) on
   the frozen feature pipeline. DISCOVERY threshold: hgb score in the
   top 1% of ALL funnel-surviving (peptide, allele) scores in this run
   (cohort-internal percentile; threshold value computed once, locked in
   the results JSON before any named-candidate selection).
7. HLA panel (locked): HLA-A*02:01, A*02:03, A*03:01, A*11:01, A*31:01,
   A*68:02, B*07:02, B*15:01 (the 8 most assay-covered alleles in the
   training data, same as the 09b G1 panel).

## Discovery criterion (from prereg, executable form)
A NAMED discovery requires ALL of: (1) zero exact matches in the locked
novelty catalogs at the disclosed lock timestamp; (2) hgb score above the
locked top-1% threshold; (3) carrier-sample TPM >= 1 AND carrier TPM >
GTEx multi-tissue median TPM (both reported); (4) MHCflurry percentile
<= 2.0 for the reported allele(s) (v1 score reported alongside);
(5) >= 1 independent evidence layer: (a) PRIDE/CPTAC mass-spec record of
the exact peptide, or (b) PubMed record (NCBI eutils) naming the exact
mutation or peptide. Zero passes = documented negative, preserved.

## Compute discipline
Streaming per-aliquot aggregation (2GB box); checkpoints committed+pushed
every ~30 min of compute; resume-safe scripts; raw manifests + thresholds
+ exclusion counts committed under results/discovery/.
