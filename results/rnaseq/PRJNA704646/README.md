# PRJNA704646 results

HelixForge completed the paired-end RNA-seq `import` path for
[PRJNA704646](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA704646): 24 technical
runs were assigned to eight biological sample pools, then processed through
QC, Salmon, tximport and MultiQC against SM_V10/WBPS19. The run used the
`v1.1.0` scientific configuration plus the merged Import metadata fix
([HelixForge PR #97](https://github.com/GMiguelAlves/HelixForge/pull/97),
commit `9fb09ff`). All 48 FASTQs passed the integrity gate. The completed
Nextflow trace contains 171 successful tasks: 169 reused from cache and two
Import tasks executed after the fix.

The eight samples produced count, TPM and effective-length matrices for 9,914
genes. Three adult pools are flagged `REVIEW` by the predefined operational
Salmon mapping threshold of 60% (53.74–58.07%); all eight passed trimming
retention, and no sample was excluded automatically. See the
[execution review](reports/PRJNA704646_execution_report.md) and
[per-sample QC](qc/qc_sample_summary.tsv).

This study is **descriptive only**. Each resistance × stage × sex pool has one
biological unit, despite its three technical runs. Treating those runs as
independent replicates would be pseudoreplication; no DESeq2 model or
differential-expression claim was made.

The public package includes standardized QC, the three matrices, the native
Import API manifest, reference provenance, a MultiQC copy with private work
paths removed, and [publication validation](manifests/publication_validation.json).
Import mode does not emit a terminal RNA-seq Integration API manifest.
Complete re-entry evidence, FASTQ checksums and command logs are preserved
privately; FASTQs, workdirs, per-read FastQC and complete Salmon directories
are excluded from Git.

`TECHNICAL_EXECUTION = PASS`

`QC = REVIEW_3_OF_8`

`SCIENTIFIC_CLASSIFICATION = DESCRIPTIVE_ONLY`
