# PRJNA551487 results — excluded from atlas

HelixForge `v1.1.0` completed the paired-end RNA-seq workflow for
[PRJNA551487](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA551487): five
DMSO and five FPL biological libraries. All 20 FASTQs passed integrity
checks, and 117/117 Nextflow tasks completed against SM_V10/WBPS19. This
compact public package contains standardized QC, gene-level counts/TPM/
lengths for 9,914 genes, sanitized MultiQC, replicate diagnostics and
provenance.

**The differential-expression outputs are withheld, not accepted biological
results.** Salmon mapping was 12.41–15.02% in all ten libraries, below the
predefined 60% review threshold. The DESeq2 process executed, but its output
is retained privately and is not published or promoted to the atlas. A
prespecified, whole-file random sample of 3,000 pairs per library found
86.60–92.87% of sampled pairs compatible with a 50-kb repetitive genomic
locus under a permissive local-alignment screen. This is **not** a molecular
identification of rRNA or proof of a library-preparation mechanism. The
published study also collected controls and treated worms from different
anatomical sites, confounding the treatment comparison. No sample was
automatically excluded or relabeled. See the [execution review](reports/PRJNA551487_execution_report.md)
and [per-library diagnostic summary](reports/fullfile_locus_qc_summary.tsv).

The [terminal manifest](manifests/rnaseq_run_manifest.json) records the full
private run and was validated against its original result tree. It is
provenance, not a claim that every private artifact is distributed in this
compact public package. The [publication record](manifests/publication_validation.json)
identifies the withheld DE outputs. FASTQs, complete Salmon directories,
individual FastQC, workdirs and scheduler logs remain outside Git.

`TECHNICAL_EXECUTION = PASS`

`SAMPLE_QC = REVIEW_10_OF_10`

`DE_PROMOTION = WITHHELD_AFTER_REVIEW`

`ATLAS_INCLUSION = EXCLUDE_FROM_ATLAS`
