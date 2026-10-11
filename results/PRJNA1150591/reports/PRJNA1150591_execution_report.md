# PRJNA1150591 — RNA-seq execution review

## Scope

[PRJNA1150591](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1150591)
contributes 15 paired-end libraries in three recorded conditions (five
biological libraries each): 7-week adult control, 7-week adult meclonazepam,
and 4-week juvenile meclonazepam. HelixForge `v1.1.0`, Nextflow `25.10.7`,
Salmon `1.10.3` and the SM_V10/WBPS19 reference were used. The native
QC → Salmon → Import → DESeq2 path finished 167/167 tasks. All 30 FASTQ
files passed the input integrity checks.

## Technical QC

| Measure | Observed |
| --- | ---: |
| Samples classified technical `PASS` | 15/15 |
| Minimum post-trim read retention | 99.67% |
| Minimum Salmon mapping | 67.79% |
| MultiQC | `PASS` |
| Gene-level import | 9,914 genes × 15 samples |
| Prespecified adult contrast completed | 1/1 |

The review triggers were read retention below 80% or Salmon mapping below
60%. Neither fired. The [QC summary](../qc/qc_summary.json) and
[sample table](../qc/qc_sample_summary.tsv) contain the underlying values.

## Statistical scope

The frozen [DE specification](../../../config/PRJNA1150591/de_spec.json)
uses `~ condition`, Wald testing, `alpha = 0.05` and a log2-fold-change
threshold of 1. All 15 libraries, including the five juveniles, were kept in
the model. The only extracted contrast is 7-week adult meclonazepam versus
7-week adult control. Retaining juveniles does **not** turn this into an age
contrast: no juvenile-versus-adult effect is reported. This choice preserves
the whole study in the fitted model and is documented explicitly rather than
silently omitting a condition.

The [contrast table](../differential_expression/contrasts/DEG_condition__adult_7wk_mclz_vs_adult_7wk_control.tsv)
contains the complete adult comparison. Of 9,458 tested genes, 936 were
classified significant by the frozen output rule. The
[DE summary](../differential_expression/deg_summary.tsv) and
[DE manifest](../differential_expression/de_manifest.json) record the model
and counts. Biological interpretation of individual genes and the juvenile
substructure seen in exploratory ordination remain pending; no samples were
relabelled or excluded.

## Provenance and limits

The terminal RNA-seq manifest passed schema, semantic and original-filesystem
validation. Public matrices and figures are checksummed in `SHA256SUMS`;
complete re-entry evidence is preserved separately in private project
storage. The public package intentionally omits raw reads, full Salmon
quantifications, RDS objects and runtime internals. This is a completed
technical execution with biological review still open.
