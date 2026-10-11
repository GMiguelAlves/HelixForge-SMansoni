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
filesystem validation. This original 20-sample result is technically complete,
but its contrast tables are retained **for audit only** following the separate
non-adult refit described below.

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

No automatic relabelling or sample exclusion was performed in the original
workflow. The adult female-versus-male contrast is `NOT_EVALUABLE`. Original
20-sample contrast tables are provided for audit, not promoted as accepted
atlas evidence. The next adult-specific review should examine independent
sex/stage markers and original sample provenance without redefining labels
from the PCA alone.

## Non-adult refit and scientific decision

A separate Slurm DESeq2 1.42.1 refit used the validated import counts and
sample table, retaining 16 non-adult samples. All four adult libraries were
omitted from this refit only; original scientific settings and outputs were
preserved. Four female-versus-male stage contrasts were computed over 9,501
filtered genes. The eight non-adult replicate pairs had VST correlations
0.945–0.994. No genes significant in both fits reversed effect direction.

| Stage | Original significant | Refit significant | Retained original | Set Jaccard |
| --- | ---: | ---: | ---: | ---: |
| Cercaria | 3,692 | 3,878 | 99.57% | 0.944 |
| Schistosomulum 1 | 147 | 342 | 91.84% | 0.381 |
| Schistosomulum 2 | 795 | 1,115 | 97.74% | 0.686 |
| Schistosomulum 3 | 207 | 454 | 94.69% | 0.422 |

The significant-gene sets, particularly stages 1 and 3, are sensitive to the
model scope. Therefore the four new non-adult contrasts are
`PASS_WITH_LIMITATIONS`, not interchangeable with the original 20-sample
tables. See [the refit contract](../reanalysis/nonadult_16/README.md) and
[comparison metrics](../reanalysis/nonadult_16/comparison_to_original.tsv).
The atlas has not yet been rewired to this supplement and must not ingest the
original tables as accepted inferential evidence.

```ini
WORKFLOW_EXECUTION = PASS
FASTQ_INTEGRITY = PASS
TECHNICAL_QC = PASS
FRAGMENT_PRIOR_SENSITIVITY = PASS
ADULT_REPLICATE_CONCORDANCE = REVIEW
ORIGINAL_20_SAMPLE_FIT = AUDIT_ONLY
NONADULT_16_SAMPLE_REFIT = PASS_WITH_LIMITATIONS
ADULT_FEMALE_VS_MALE = NOT_EVALUABLE
ATLAS_INCLUSION = PENDING_NONADULT_WIRING
```
