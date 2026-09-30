# Results policy

Every study exposes the same compact public interface under
`results/<study>/`:

- `qc/`: `qc_sample_summary.tsv`, `qc_summary.json`, and optionally MultiQC;
- `expression/`: gene-level count, abundance and effective-length matrices;
- `differential_expression/`: model summaries, accepted contrasts and figures
  when differential expression is scientifically supported;
- `manifests/`: the terminal RNA-seq manifest and its validation record;
- `reports/`: execution, validation, performance and interactive reports;
- `README.md`: study scope, status and limitations.

`gene_report/` is optional because it is a presentation product. A missing
`differential_expression/` payload is valid only when the study README records
that inference was not performed by design.

The public package is not the execution workdir and is not the private
re-entry archive. Raw FASTQs, per-read FastQC reports, Salmon auxiliary trees,
complete `quant.sf` collections, duplicated Integration API payloads, RDS
objects, scheduler logs, cache, references and indexes remain in verified
persistent storage or audit archives. They must not be copied into Git merely
to make a public terminal manifest self-contained.
This exclusion applies equally to every registered study; directories such as
`compatibility/quants/` are not part of the public results contract, including
for historical projects.

The stable paths consumed by reviewers and the UI are therefore the compact
directories above. Native Nextflow output trees may have a different internal
layout, but publication must map their accepted products onto this interface.
Study-specific finalizers must not introduce a new public directory layout.

The integrated descriptive RNA-seq atlas is available under `atlas/rnaseq/`.
It contains a navigable HTML report, SVG figures, machine-readable tables,
browser data, and a checksum manifest. It does not include raw reads, private
runtime paths, or a pooled cross-study inferential model.
