# E-ERAD-478 RNA-seq execution report

## Classification

`PASS`

The complete paired-end RNA-seq workflow finished successfully under
HelixForge v1.0.2. All scientific outputs were reviewed before cleanup and the
accepted result and audit packages were preserved in persistent storage.

## Scope

- 120 ENA runs;
- 60 biological samples;
- two technical runs aggregated per biological sample;
- WBPS19 / SM_V10 reference bundle;
- Salmon quantification followed by tximport;
- DESeq2 design `~ batch + condition`;
- 32 prespecified contrasts.

## Validation summary

- workflow: 950 of 950 processes completed, zero failed;
- quantification: 60 of 60 biological samples;
- QC: 60 `PASS`, zero `REVIEW`, zero `FAIL`;
- Salmon mapping: 75.03% to 91.84%, median 85.74%;
- expression matrices: 9,914 genes and 60 sample columns;
- differential expression: all 32 frozen contrasts completed;
- terminal manifest: validated;
- MultiQC and candidate-gene report: completed;
- cleanup: completed only after persistence and audit verification.

## Public and private products

The Git package contains the compact QC interface, gene-level matrices,
differential-expression products, reports, figures and manifests. Complete
Salmon output directories, duplicated Integration API artifacts, RDS objects,
per-read FastQC reports and the full candidate-gene report remain in the
verified private storage/audit package.

No scientific result was modified during public packaging. The package only
maps accepted outputs to the repository's stable results contract and removes
reproducible duplication from Git.
