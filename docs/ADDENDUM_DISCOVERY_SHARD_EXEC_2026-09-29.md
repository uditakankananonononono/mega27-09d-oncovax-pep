# 09d prescreen compute-sharding addendum - 2026-09-29

Locked before remote worker scores: user requested faster parallel execution.
This is an execution-only change, not a new cohort, allele panel, outcome,
model, threshold or candidate selection. The original process keeps its
contiguous on-box prefix. Twelve independent GitHub Actions workers score
non-overlapping chunk ranges 2581-3284. Chunk inputs are sorted unique
peptides from the frozen self-filtered records; SHA256 of decompressed lines
must equal `cf9b52a3a997e1c07a76fa5168375aad4638a5317e6d054dfa6deaf76bf82064`
with 8,212,079 peptides. Actual count is 3,285 chunks (2,500 peptides per
full chunk), correcting earlier progress denominator 3,278. No chunk is
accepted twice. Workers retain MHCflurry 2.2.1 and the same model download,
8-allele order, 20,000-pair chunk boundaries, percentile gate <=2.0, and
serialized output schema. Runtime version/device variance is checked against
an existing box chunk before merger. Workers push disjoint branches, not main;
main merge requires range completeness, hashes, row count and pair order
checks. The box must stop before scoring the worker-assigned suffix to avoid
racing and double counting. All raw negative and exclusion counts remain.

Box boundary is enforced by PRESCR_BOX_END_CHUNK=2581, returning before
any cross-shard chunk. This is an execution guard only; no QC is emitted
until all worker chunks are merged and the standard prescr phase is resumed.

At the handoff, the box completed chunk 2580 while being stopped; its full
20,000-pair output was validated against the frozen peptide boundary and
committed at e642583. Remote range therefore starts at 2581, not 2580.
The first workflow dispatched at cebf417 starts 2580 and must NOT be merged
for chunk 2580; prefer the corrected workflow run and verify no duplicate.

Worker branch publication failed in the first run after scoring succeeded;
its computed files were lost with ephemeral runners. The recovery workflow
uploads each range as a named artifact before the job ends. Branch publication
is removed. A run with both artifact upload and branch push may finish with
an error even when artifact exists; the artifact's bytes and peptide-pair
alignment, not the job success badge, determine whether chunks are usable.
