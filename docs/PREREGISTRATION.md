# PRE-REGISTRATION DRAFT - mega27-09d (working title: OncoVax-Pep)
Status: LOCKED 2026-09-26 17:05 IST, before any outcome data was scored.
Source: ChatGPT ideation round 1 (2026-09-26, conversation 6ab7a277), independently verified for dataset/comparator existence; user rules 1-7 apply.

## Question
Can a multimodal (sequence + allele + expression + processing) computational pipeline identify tumor neoantigen peptides that are experimentally immunogenic, ranking them better than binding-affinity-only pipelines, and surface NAMED, previously uncatalogued candidates that survive independent evidence layers?

## Winner archetype (rule 7)
Identify + verify (Kulviwat claudin-5 archetype, verified at isef.net/project/bmed028-...): candidate discovery from public data, then layered independent verification (expression, HLA presentation, proteomics, pathway, structure).

## Data (public, verified to exist 2026-09-26)
- IEDB T-cell assay export (immunogenic / non-immunogenic labels, peptide + HLA allele): iedb.org database export v3.
- TCGA somatic mutations: GDC open-access MAFs (docs.gdc.cancer.gov cohort MAF).
- TCGA expression: UCSC Xena (xena.ucsc.edu).
- Proteomic validation (optional layer): CPTAC/PRIDE.

## Locked splits
Allele-group + peptide-cluster holdout (no peptide cluster or allele family shared between train/test). Frozen before training.

## Benchmark gate (beat, matched evaluation)
- Comparator A: NetMHCpan-4.1 binding-rank baseline on identical immunogenicity test set (published tool, run via free IEDB tools or published metrics on identical partition).
- Comparator B: a published immunogenicity predictor (DeepImmuno, PMC7781330; NeoTImmuML, PMC12585993) on identical test partition.
- Metrics (locked): AUPRC, Recall@Top-100, NDCG (NOT AUROC alone). BEAT = exceeds both comparators on the primary metric (AUPRC) on the identical frozen test set; fold-level paired test vs our own ablations.

## Discovery criterion (locked, falsifiable)
A candidate counts as a NEW DISCOVERY only if ALL hold: (1) not present in IEDB/neoantigen catalogs at lock time; (2) top-1% model score; (3) tumor-specific expression support in TCGA cohorts; (4) predicted HLA presentation; (5) at least one independent evidence layer (proteomics or literature). If zero candidates pass: documented negative, pivot per rule 6.

## Honest-negative protocol
All failed gates preserved in results/ with the same voice; no re-fishing after outcomes are seen; amendments logged as addenda, never edits.
