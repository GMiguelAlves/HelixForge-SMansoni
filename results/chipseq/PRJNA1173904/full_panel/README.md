# Full eight-mark panel

Accepted with warnings: 64 IP libraries, eight marks, four gonad/mating conditions,
two labelled replicates per combination. See the [study overview](../README.md)
and [processing policy](../../../../docs/chipseq/prjna1173904_analysis.md).

- [Per-library calls](peaks/), BED6 gzip, zero-based half-open, no consensus recutting.
- [HOMER annotations](annotation/), including library identity, control policy,
  QC/review warnings and gene-assignment ambiguity fields.
- [Samples](metadata/sample_table.tsv), [QC](qc/qc_sample_summary.tsv),
  [model metrics](metrics/model_summary.tsv), [annotation metrics](metrics/annotation_summary.tsv),
  [replicate overlap](metrics/replicate_region_overlap.tsv), and [fit figures](figures/).
- [Parameters](manifests/parameters.json) and [versions](manifests/software_versions.tsv).

Counts summed across libraries are not unique loci. Testis results are a
no-Input alternative, not a controlled reproduction. Acceptance with warnings
does not establish every called interval as true occupancy or every nearest-gene
association as regulatory. The native ChromstaR BED score is not a MACS p/q.
Gene ambiguity flags assess nearest-TSS distance ties, not all biological ambiguity.

One compact annotated table per library includes the nearest-gene assignment; a second duplicate association table is not distributed. Join record_id to the sample/QC tables for full library metadata. Projection preserves every interval, gene assignment, control policy, warning and ambiguity field.
