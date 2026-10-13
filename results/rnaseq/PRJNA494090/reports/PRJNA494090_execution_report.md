# PRJNA494090 — execution and QC review

## Completed path

The paired-end HelixForge `v1.1.0` run completed 259/259 tasks for 24
biological libraries (P6/P7/P8 × selected/control × four replicates). All
48 FASTQs passed size, checksum, gzip and pairing validation. Salmon 1.10.3,
tximport, MultiQC and the three frozen DESeq2 contrasts executed with the
SM_V10/WBPS19 reference. The count, TPM and effective-length matrices contain
9,914 genes × 24 libraries; their checksums match the native
[Import manifest](../expression/import_manifest.json).

## Independent QC gate

| Measure | Observed |
| --- | ---: |
| Libraries classified `PASS` | 14/24 |
| Libraries classified `REVIEW` | 10/24 |
| Post-trim retention | 97.72–99.48% |
| Salmon mapping | 49.28–84.35% |
| Salmon mapping mean | 64.53% |
| Minimum within-group median replicate correlation | 0.9718 |
| MultiQC | `PASS` |

The predefined review triggers are post-trim retention below 80% and Salmon
mapping below 60%. Only the mapping trigger fired. Group means were 53.35%
(P6 control), 63.06% (P6 selected), 62.42% (P7 control), 79.13% (P7
selected), 55.32% (P8 control) and 73.92% (P8 selected). This asymmetry
tracks the experimental comparison and could bias interpretation. Replicate
correlation and the [exploratory PCA](replicate_pca.png) show consistency of
the mapped counts, but cannot rule out systematic differences in the
unmapped fraction. The [per-sample QC](../qc/qc_sample_summary.tsv) and
[replicate summaries](replicate_group_summary.tsv) provide the underlying
values.

## Scientific decision

The libraries are mixed-sex, and the study may contain differences in sex
composition across selected and control material. The QC review and sample
composition were considered independently of the number of DE genes.
Therefore the DESeq2 outputs are **not promoted as biological results** and
remain in private audit storage. No contrast table is included here, no
sample has been excluded, and no unmeasured sex covariate was added.

Further work should characterize the unmapped fraction and sample-level sex
composition before interpreting resistance-associated effects. The
[terminal manifest](../manifests/rnaseq_run_manifest.json) and its
[validation](../manifests/terminal_manifest_validation.json) refer to the
complete private execution; the compact public package intentionally omits
DE artifacts while promotion is on `HOLD`.
