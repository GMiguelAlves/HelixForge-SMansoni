# PRJEB14695 results

Accepted compact outputs from the full native RNA-seq execution of 138
paired-end technical runs representing 23 biological samples.

- `expression/`: gene-level counts, TPM, effective lengths, import metadata,
  versions, and R session information.
- `differential_expression/`: the frozen DESeq2 model, six contrasts,
  normalized counts, aggregate tables, and scientific figures.
- `manifests/` and `rnaseq/`: portable terminal manifests and integration-ready
  artifacts.
- `qc/`: compact QC and Salmon mapping summary.
- `reports/`: MultiQC, Nextflow HTML reports, validation, performance,
  sanitization evidence, and the reviewed execution report.

FastQC HTML files and heavy intermediates are intentionally excluded. The
public results contain no private server paths, usernames, job identifiers, or
raw sequencing data.
