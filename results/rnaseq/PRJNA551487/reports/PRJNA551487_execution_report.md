# PRJNA551487 — execution and QC review

## Completed path

The paired-end HelixForge `v1.1.0` run completed 117/117 tasks for ten
biological libraries (five DMSO, five FPL). All 20 FASTQs passed size,
checksum, gzip and pairing validation. Salmon 1.10.3, tximport, MultiQC and
the frozen DESeq2 contrast executed with the SM_V10/WBPS19 reference. The
count, TPM and effective-length matrices contain 9,914 genes × ten
libraries; their checksums match the native
[Import manifest](../expression/import_manifest.json).

## Independent QC gate

| Measure | Observed |
| --- | ---: |
| Libraries classified `PASS` | 0/10 |
| Libraries classified `REVIEW` | 10/10 |
| Post-trim retention | 99.75–99.81% |
| Salmon mapping | 12.41–15.02% |
| Salmon mapping mean | 13.96% |
| Mapped/processed fragments | 23,514,623/168,740,503 |
| MultiQC | `PASS` |

The predefined review triggers are post-trim retention below 80% and Salmon
mapping below 60%. Only the mapping trigger fired. The mapped fraction was
12.41–13.81% in DMSO and 14.70–15.02% in FPL. The
[exploratory PCA](replicate_pca.png) and
[replicate correlations](replicate_correlations.tsv) describe the mapped
counts, but cannot rule out bias in the much larger unmapped fraction.

A small independent BLAST-based diagnostic found frequent strong hits to a
repetitive rDNA locus in reads sampled from the beginning of the FASTQs.
The sampling did **not** cover each complete library, so those hit fractions
must not be presented as whole-library rDNA prevalence. The specific
library-preparation cause remains unproven. A whole-file, pre-specified
sampling design would be needed before making a quantitative attribution.

## Scientific decision

The QC decision was recorded independently of differential-expression
results. The DESeq2 outputs are **not promoted as biological results** and
remain in private audit storage. No contrast table is included here, no
sample has been excluded, and no threshold was changed. The original study
also collected DMSO controls and FPL-treated worms from different anatomical
sites, so treatment and collection site cannot be separated in this design.

The [terminal manifest](../manifests/rnaseq_run_manifest.json) and its
[validation](../manifests/terminal_manifest_validation.json) refer to the
complete private execution; the compact public package intentionally omits
DE artifacts while promotion is on `HOLD`.
