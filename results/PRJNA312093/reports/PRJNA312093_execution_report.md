# PRJNA312093 — single-end RNA-seq execution review

## Scope and result

Picard et al. [PRJNA312093](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA312093)
contributes 20 single-end libraries: female and male samples at five recorded
stages, with two biological replicates per group. The frozen design estimates
five female-versus-male contrasts within stage using `~ condition`. HelixForge
`v1.1.0` (commit `062001ceb46d8a25c2415823b97d0e310658f8c0`), Nextflow
`25.10.7` and Java `21` executed the native QC → Salmon → Import → DESeq2 path
with the SM_V10/WBPS19 reference and a validated prebuilt Salmon index.

The run completed 161/161 scientific tasks. All 20 input FASTQs passed MD5
and gzip checks. The terminal manifest passed schema, semantic and original
filesystem validation. The result is technically complete, but **biological
promotion remains on hold**.

## Prespecified technical QC

| Measure | Observed result |
| --- | ---: |
| Samples classified technical `PASS` | 20/20 |
| Minimum post-trim read retention | 88.78% |
| Minimum Salmon mapping | 67.59% |
| MultiQC | `PASS` |
| Gene-level import | 9,914 genes × 20 samples |
| DESeq2 contrasts completed | 5/5 |

The prespecified review triggers were retention below 80% or Salmon mapping
below 60%. Neither trigger fired. In the single-end QC summary, historical
`*_both_mates` column names mean R1 counts only.

## Single-end fragment-prior sensitivity

The frozen primary Salmon prior was `--fldMean 200 --fldSD 80`. A 250/80
alternative was selected before examining DE results and executed separately
without replacing the primary output. The primary `tximport` counts reproduced
the workflow output exactly. Among 9,526 genes common to both statistical
models, the five contrasts retained 95.17–99.65% of the primary significant
genes; set Jaccard ranged from 0.902 to 0.994. No gene significant under both
priors reversed direction. These are stability checks, not a method for
choosing the prior with more significant genes. Detailed alternative results
and checksums are retained in the private audit package.

## Biological replicate review

VST correlation for the two samples labelled adult male was 0.800; the other
nine replicate pairs ranged from 0.943 to 0.994. Each adult-male sample was
more similar to a cercaria-labelled sample than to its designated replicate.
The published run accessions, frozen mapping, downloaded file checksums and
Salmon input commands were checked; they do not show an internal file swap.
They cannot establish the original biological identity of the libraries.

No automatic relabelling or sample exclusion was performed. The adult
female-versus-male contrast is especially sensitive to this unresolved
discordance. Other contrast tables are provided for audit, not promoted as
accepted atlas evidence. The next review should examine independent
sex/stage markers and original sample provenance without redefining labels
from the PCA alone.

```ini
WORKFLOW_EXECUTION = PASS
FASTQ_INTEGRITY = PASS
TECHNICAL_QC = PASS
FRAGMENT_PRIOR_SENSITIVITY = PASS
BIOLOGICAL_REPLICATE_CONCORDANCE = REVIEW
STUDY_PROMOTION = HOLD
ATLAS_INCLUSION = HOLD
```
