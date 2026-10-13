# PRJNA551487 results — QC review hold

HelixForge `v1.1.0` completed the paired-end RNA-seq workflow for
[PRJNA551487](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA551487): five
DMSO and five FPL biological libraries. All 20 FASTQs passed integrity
checks, and 117/117 Nextflow tasks completed against SM_V10/WBPS19. This
compact public package contains standardized QC, gene-level counts/TPM/
lengths for 9,914 genes, sanitized MultiQC, replicate diagnostics and
provenance.

**The differential-expression outputs are on scientific `HOLD`, not accepted
results.** Salmon mapping was 12.41–15.02% in all ten libraries, below the
predefined 60% review threshold. The DESeq2 process executed, but its output
is retained privately and is not published or promoted to the atlas. A
limited, early-FASTQ diagnostic found frequent matches to a repetitive rDNA
locus; it does not estimate whole-library composition or prove a particular
library-preparation cause. No sample was automatically excluded or
relabeled. See the [execution review](reports/PRJNA551487_execution_report.md).

The [terminal manifest](manifests/rnaseq_run_manifest.json) records the full
private run and was validated against its original result tree. It is
provenance, not a claim that every private artifact is distributed in this
compact public package. The [publication record](manifests/publication_validation.json)
identifies the withheld DE outputs. FASTQs, complete Salmon directories,
individual FastQC, workdirs and scheduler logs remain outside Git.

`TECHNICAL_EXECUTION = PASS`

`SAMPLE_QC = REVIEW_10_OF_10`

`DE_PROMOTION = HOLD`

`ATLAS_INCLUSION = HOLD`
