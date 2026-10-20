# PRJEB32839 results

Accepted compact records from the full native RNA-seq execution of 150
paired-end technical runs representing 75 biological samples.

- `reports/`: reviewed execution report and compact report figures.
- `qc/`: per-sample QC classification and aggregate trimming/mapping metrics.
- `differential_expression/`: the 12 frozen contrast summaries.
- `gene_report/`: compact candidate-gene tables and overview figures.
- `manifests/`: portable terminal manifest and its validation record.

The complete 1.19 GB result tree, including 750 FastQC reports, MultiQC,
Nextflow HTML reports, expression matrices, complete DE tables, and the full
candidate-gene HTML report with 348 figures, is preserved in verified private
persistent storage. Raw FASTQs, work directories, references, indexes, cache,
and private cluster logs are intentionally excluded from Git.

`PRJEB32839_RNASEQ_ANALYSIS = PASS_WITH_LIMITATIONS`

The study-specific scratch tree was removed only after persistent-copy and
audit gates passed. Cleanup recovered 1.208 TB, left no heavy residuals, and
did not touch shared reference/runtime resources or other studies.
