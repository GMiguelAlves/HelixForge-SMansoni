# PRJEB3190 results

The complete native paired-end RNA-seq execution is classified `PASS`. The
workflow completed 227/227 processes without failures and reused the certified
WBPS19/SM_V10 Salmon index.

- 20 paired-end biological samples spanning seven schistosomulum time points
  were processed.
- Salmon produced 20 sample quantifications.
- tximport produced gene-level count, TPM and effective-length matrices with
  9,914 genes and 20 sample columns.
- DESeq2 used the frozen `~ condition` design and produced nine prespecified
  contrasts: six consecutive time-point comparisons and three comparisons to
  the 0 h baseline.
- The optional sequencing-batch labels remain exploratory because their ENA
  aliases were not documented as technical batches; they were not included in
  the inferential model.
- MultiQC, the candidate-gene report, Nextflow reports and the terminal RNA-seq
  manifest completed successfully.

Directory guide:

- `qc/`: consolidated MultiQC report.
- `expression/`: gene-level matrices, import statistics and R provenance.
- `differential_expression/`: normalized counts, contrast tables and figures.
- `reports/`: Nextflow execution reports and performance trace.
- `rnaseq/`: portable terminal manifest and Integration API artifacts.
- `manifests/`: terminal validation, reference and quantification manifests.
- `compatibility/quants/`: exact Salmon `quant.sf` files preserved for re-entry.

The complete candidate-gene report was validated and is retained in the
project's persistent private result storage. It is intentionally excluded from
Git to keep the public repository compact.

The authoritative terminal manifest is
[`rnaseq/rnaseq_run_manifest.json`](rnaseq/rnaseq_run_manifest.json).

`PRJEB3190_RNASEQ_ANALYSIS = PASS`
