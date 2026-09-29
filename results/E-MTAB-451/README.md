# E-MTAB-451 results

The complete native paired-end RNA-seq execution is classified `PASS`. The
workflow completed 136/136 processes without failures and reused the certified
WBPS19/SM_V10 Salmon index.

- 12 paired-end runs representing 11 biological samples were processed.
- The two technical runs assigned to `schistosomulum_3h_r1` were consolidated
  before quantification.
- Salmon produced 11 sample quantifications; mapping rates ranged from 62.26%
  to 78.87%, with a median of 75.18%.
- tximport produced gene-level count, TPM and effective-length matrices with
  9,914 genes and 11 sample columns.
- DESeq2 used `~ condition` and produced the three prespecified contrasts.
- `adult_mixed` and `cercarial_tail` each contain one sample and remain
  descriptive; they were not used in inferential contrasts.
- MultiQC, the candidate-gene report, Nextflow reports and the terminal RNA-seq
  manifest completed successfully.

Directory guide:

- `qc/`: consolidated MultiQC report.
- `expression/`: gene-level matrices, import statistics and R provenance.
- `differential_expression/`: normalized counts, contrast tables and figures.
- `reports/`: Nextflow reports, candidate-gene HTML bundle and trace.
- `rnaseq/`: portable terminal manifest and Integration API artifacts.
- `manifests/`: terminal validation, reference and quantification manifests.
- `compatibility/quants/`: exact Salmon `quant.sf` files preserved for re-entry.

The authoritative terminal manifest is
[`rnaseq/rnaseq_run_manifest.json`](rnaseq/rnaseq_run_manifest.json).

`E_MTAB_451_RNASEQ_ANALYSIS = PASS`

