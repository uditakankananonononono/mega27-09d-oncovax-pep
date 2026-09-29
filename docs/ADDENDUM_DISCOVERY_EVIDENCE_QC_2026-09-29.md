# 09d candidate evidence QC, 2026-09-29

This addendum records a post-threshold source check. It does not change the locked HGB threshold, re-rank candidates, or establish an immunogenic discovery.

## GVFGGLWGV source and exact-self check

The frozen reviewed-human-proteome index contains CLCN3 canonical sequence (818 aa), with A at one-based position 405. The stored mutation record `TCGA-EE-A3AD-06A`, `CLCN3 p.A405V` produces `GVFGGLWGV` at one-based positions 396-404. The corresponding WT window is `GVFGGLWGA`. The mutant peptide is present in `neopeptides.jsonl.gz` and `funnel_self.jsonl.gz`. Direct exact substring search found no `GVFGGLWGV` among the 20,213 sequences of the stored index or within the raw stored reviewed-human FASTA. This validates the exclusion only against the declared frozen reference, not against every isoform, variant or population proteome.

PeptideAtlas ProteoMapper's mapping endpoint (`https://peptideatlas.org/api/promast/v1/map?proteome=Hs&peptide=GVFGGLWGV`) reports one mapping to P04637 at 273, even with `fuzzy=0`. This result is inconsistent with the current UniProt P04637 canonical TP53 sequence (`https://rest.uniprot.org/uniprotkb/P04637.json`), which has no exact `GVFGGLWGV`; its positions 273 onward do not contain that sequence. The discrepancy's cause is unresolved, including the PeptideAtlas database version and mapping semantics. The returned mapping must not be treated as an exact reference-self hit, nor as measured peptide presentation. No spectrum-level identification has been established by this check. Until the discrepancy is resolved, do not claim exhaustive self exclusion.

The current lead has HGB score 0.9712772070663221, MHCflurry affinity 13.20263682606981 nM, percentile 0.02475 and source-sample tumor TPM 273.6284. These are in-silico and transcript evidence, not immunopeptidomics confirmation. GTEx v8 CLCN3 median is 16.8262 TPM across 54 tissues. Exact gene/variant PubMed and PRIDE project-keyword searches did not return hits, but are not peptide-level absence tests.

## Next independent evidence gate

A candidate can only be described as experimentally presented if an independent peptide-level identification with a traceable dataset, spectrum/PSM, sample and sequence is found and reviewed. A peptide-to-protein mapping or project-keyword result is not that evidence. If no identification is found, label the candidate computational and list the unresolved evidence explicitly.

## Independent mutation-catalog cross-check

CAN-IMMUNE's live CLCN3 mutation endpoint (`https://canelib.erc.monash.edu/api/gene_mutations?gene=CLCN3&search%5Bvalue%5D=A405V&length=10`, checked September 29) returns two transcript records for `p.A405V`, `c.1214C>T`, `TCGA-EE-A3AD-06`, Skin, source `cosmic_TS`. Both list mutant local peptide `FILLGVFGGLWGVFFIRANIAWCRR` and wild peptide `FILLGVFGGLWGAFFIRANIAWCRR`. Its sample label lacks the final `A` seen in the analysis input (`TCGA-EE-A3AD-06A`); this is a sample-alias difference to reconcile rather than silently conflate. This corroborates the variant and mutant local sequence independently of our generation script, but is not a peptide-spectrum match or HLA presentation experiment.

## HLA Ligand Atlas benign reference

The official 2020.12 archive (`https://hla-ligand-atlas.org/rel/hla_2020.12.zip`; SHA256 `e9622ddaab5cb18592e9378b5df917b73afcbe6e08a5151be99b18aa59a232d1`) contains 223,246 unique sequence rows in its aggregated table. Exact sequence join against all 7,007 unique peptides in the locked 09d shortlist found zero matches, including the top 40 ranked unique peptides. Reproduce with `python3 scripts/discovery_hla_atlas_check.py /path/to/hla_2020.12.zip`; machine-readable result is `results/discovery/hla_atlas_2020_12_exact_check.json`. The release's table semantics and license are explained at `https://hla-ligand-atlas.org/data`. A no-hit in this benign-tissue release is neither tumor presentation evidence nor proof that the peptide was never measured in healthy tissue.

## dbPepNeo 2019 historical neoantigen catalog check

The HC T-cell-response and MC mass-spectrometry-plus-exome archives linked at `http://www.biostatistics.online/dbPepNeo/download.html` were downloaded and exact-matched to all 7,007 shortlisted peptide sequences. Neither has a match. Their SHA256 hashes, data-row counts, and exact intersections are recorded in `results/discovery/dbpepneo_2019_exact_check.json`; reproduce with `python3 scripts/discovery_dbpepneo_check.py /path/to/HC.zip /path/to/MC.zip`. These small 2019 catalogs cannot establish global novelty or tumor presentation. Their addition is an **exploratory post-lock cross-check**, not a retroactive change to the preregistered novelty catalogs or thresholds.
