# OncoVax-Pep data provenance (pinned)
- IEDB T-cell assay export v3: data/raw/tcell_full_v3.zip (45,075,473 bytes),
  https://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip,
  retrieved 2026-09-27; sha256 in data/SHA256SUMS.txt. Filtered to human,
  MHC class I, canonical 8-15mer peptides, Positive*/Negative labels ->
  data/processed/iedb_tcell_class1_human.csv (build: scripts_build_dataset.py;
  the 1.35GB extracted CSV is deleted after the build, re-derivable from zip).
- TCGA somatic mutations (GDC open MAFs) and TCGA expression (UCSC Xena):
  NOT yet fetched - next pipeline steps.
- Comparators (existence verified 2026-09-26): NetMHCpan-4.1 (IEDB tools /
  local MHCflurry proxy documented), DeepImmuno PMC7781330, NeoTImmuML
  PMC12585993.
