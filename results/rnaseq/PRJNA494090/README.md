# PRJNA494090 results — QC review hold

HelixForge `v1.1.0` completed the paired-end RNA-seq workflow for
[PRJNA494090](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA494090): 24
independent biological libraries across P6, P7 and P8, with selected and
control groups (four libraries per group). All 48 FASTQs passed integrity
checks, and 259/259 Nextflow tasks completed against SM_V10/WBPS19. The
public package contains standardized QC, gene-level counts/TPM/lengths for
9,914 genes, sanitized MultiQC, replicate diagnostics and provenance.

**The differential-expression outputs are on scientific `HOLD`, not accepted
results.** Ten of 24 libraries fell below the predefined 60% Salmon mapping
review threshold. Mapping also differs between control and selected groups,
while the experimental material is mixed-sex; those factors could affect the
interpretation of resistance-associated differences. The DESeq2 processes
executed, but their tables are retained privately and are not published or
promoted to the atlas pending further review. No sample was automatically
excluded or relabelled. The [execution review](reports/PRJNA494090_execution_report.md)
explains the evidence and limitations.

The [terminal manifest](manifests/rnaseq_run_manifest.json) records the full
private run and was validated against its original result tree. It is
provenance, not a claim that every private artifact is distributed in this
compact public package. The [publication record](manifests/publication_validation.json)
identifies the withheld DE outputs. FASTQs, complete Salmon directories,
individual FastQC, workdirs and scheduler logs remain outside Git.

`TECHNICAL_EXECUTION = PASS`

`SAMPLE_QC = REVIEW_10_OF_24`

`DE_PROMOTION = HOLD`

`ATLAS_INCLUSION = HOLD`
