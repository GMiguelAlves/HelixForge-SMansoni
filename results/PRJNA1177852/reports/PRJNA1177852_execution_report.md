# PRJNA1177852 — descriptive RNA-seq execution review

## Scope

[PRJNA1177852](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1177852)
contributes eight paired-end adult libraries sampled at 0, 0.25, 1, 3, 6,
12, 24 and 48 hours after meclonazepam exposure. HelixForge `v1.1.0`,
Nextflow `25.10.7`, Salmon `1.10.3` and the SM_V10/WBPS19 reference were
used. The native QC → Salmon → Import path finished 91/91 tasks; all 16
input FASTQs passed integrity checks. DESeq2 was intentionally not run.

## Technical QC

| Measure | Observed |
| --- | ---: |
| Samples classified technical `PASS` | 8/8 |
| Minimum post-trim read retention | 99.83% |
| Minimum Salmon mapping | 66.99% |
| MultiQC | `PASS` |
| Gene-level import | 9,914 genes × 8 samples |

The predefined review triggers were read retention below 80% or Salmon
mapping below 60%. Neither fired. The [QC summary](../qc/qc_summary.json)
and [sample table](../qc/qc_sample_summary.tsv) contain the underlying
values. The expression matrices are suitable for descriptive exploration,
with their checksums recorded in the native
[Import manifest](../expression/import_manifest.json).

## Scientific boundary

There is only one biological library per time point and no contemporaneous
control at each time. A time-course DESeq2 contrast would lack the replication
and control structure required by the frozen design, so no p-values, adjusted
p-values or differential-expression claims are supplied. Cross-time patterns
are descriptive and may be influenced by time and experimental factors that
cannot be separated here.

The `import` run mode does not produce a terminal Integration API manifest.
The [publication record](../manifests/publication_validation.json) verifies
the compact result package without posing as one. Original input and re-entry
evidence remain in private project storage; Git contains neither raw reads nor
full quantification payloads.
