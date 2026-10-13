# PRJNA1177852 results

HelixForge `v1.1.0` completed the paired-end RNA-seq `import` mode for all
eight time-point libraries: QC, Salmon, tximport and MultiQC. The reference
was SM_V10/WBPS19. All 16 input FASTQs passed integrity checks, all eight
samples passed the predefined technical QC thresholds, and 91/91 workflow
tasks completed.

The recorded time points are 0, 0.25, 1, 3, 6, 12, 24 and 48 hours after
meclonazepam exposure. Each time point has one biological library and lacks a
parallel control. Therefore this package is **descriptive only**: it contains
counts, TPM and length matrices but no DESeq2 model, differential-expression
table or inferential time-course claim. Apparent time trends should not be
interpreted as statistically supported treatment effects.

The compact public package includes standardized QC, matrices, the native
Import API manifest, reference provenance, sanitized MultiQC and an
[execution review](reports/PRJNA1177852_execution_report.md). The workflow's
`import` mode did **not** emit a terminal RNA-seq Integration API manifest;
the [publication validation](manifests/publication_validation.json) documents
this distinction. Complete re-entry evidence is preserved privately. FASTQs,
individual FastQC reports, complete Salmon directories, RDS objects, workdirs,
and scheduler logs are excluded from Git.

`TECHNICAL_EXECUTION = PASS`

`SCIENTIFIC_CLASSIFICATION = DESCRIPTIVE_ONLY`
