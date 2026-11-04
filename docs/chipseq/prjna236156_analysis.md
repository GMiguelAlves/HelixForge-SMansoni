# PRJNA236156: article-guided four-mark processing

## Scope and source

[Roquis et al. (2016)](https://doi.org/10.1111/mec.13555) studied chromatin
changes in *Schistosoma mansoni* following sympatric/allopatric snail-host
exposure. The accepted [result panel](../../results/chipseq/PRJNA236156/README.md)
contains 48 individual IP peaksets and 12 matched UB0 controls, spanning four
marks in cercariae and adults. Sample identities are taken from the explicit
ENA aliases and cross-checked against Appendices S2 and S3.

The goal is to obtain mark-specific intervals for later use, not to reproduce
the paper's epimutation counts, Differential Binding or EpiCSeg states.

## Fixed processing settings

| Step | Tool and settings |
|---|---|
| QC | FastQC 0.12.1; Trim Galore 0.6.10/Cutadapt 5.1, Q20 and minimum length 25; post-trim FastQC; MultiQC 1.35 |
| Alignment | Bowtie2 2.5.5, `--end-to-end --sensitive --gbar 4`; validated prebuilt SM_V10/WBPS19 index |
| Alignment selection | SAMtools 1.23.1, MAPQ ≥40, exclude unmapped reads and all records with an XS tag |
| Duplicate cap | Pyicos 2.0.7, revision `dd9b071327320395d9e30d325c885b70f7c08705`; at most four copies per chromosome/start/end/strand |
| Sampling | GNU shuf without replacement; deterministic AES-256-CTR random source; seed identifier `20261103`; coordinate sorting afterward |
| Peak calling | PeakRanger 1.18, `ranger`, p=0.0001, FDR=0.01, extension 200 bp, bandwidth 99, delta 1; matched Input |
| QC metrics | BEDTools union coverage and per-read overlap fractions for IP and matched Input |

Pyicos counts additional duplicate copies: `--duplicates 3` retains at most
four total copies. SAM-to-BED coordinate/strand handling, this cap and exact
sampling counts were checked independently. There is no complete duplicate
removal, second caller-side deduplication or automatic model fallback.

## Per-call depths

| Mark/context | Selected reads per IP and matched Input |
|---|---:|
| H3K4me3, cercariae | 16,000,000 |
| H3K4me3, cercaria Gua Fa1 | 14,730,063, all available |
| H3K4me3, cercaria Gua Fb5 | 5,925,466, all available |
| H3K4me3, adults | 5,925,466 |
| H3K27ac/H3K27me3/H4K20me1, general | 17,774,492 |
| H3K27ac, cercaria Gua Fb8 | 12,055,170, all available |
| H4K20me1, cercaria Gua Fb8 | 8,351,562, all available |

Each IP uses its own corresponding UB0 control. Equal run/depth Input subsets
are reused between marks where appropriate; they are not extra replicates.
No oversampling creates reads to fill a missing depth.

## Adaptations relative to the article

The paper's S2 specifies 16 million H3K4me3 reads for cercariae, except Gua
Fb5 at 6 million, and 7 million for adults. The remaining published targets
are 23 million for H3K27ac/H3K27me3 and 19 million for H4K20me1, with lower
Fb8 exceptions. Actual accepted depths above reflect available reads after
the current reference, QC and selection. Input-depth cells are blank in S2;
matching the IP and Input depths is an explicit analysis policy, not a
verified historical parameter.

SM_V10/WBPS19, the standard QC tools, current Bowtie2/SAMtools, Pyicos 2.0.7,
PeakRanger 1.18 rather than 1.16, fixed MAPQ40 within the reported 40–42 range,
and the deterministic sampling implementation are recorded adaptations.
The historical Perl sampling script/seed is not supplied. This analysis
therefore does not claim exact executable or numerical reproduction.

## Review and accepted use

The [paper comparison](../../results/chipseq/PRJNA236156/metrics/paper_comparison.tsv)
retains all 48 matches and the original S3 source cell. H3K4me3 is close in
region count after depth adjustment. Other marks, especially the two cercaria
Gua Fb8 libraries flagged in S2, remain heterogeneous. More peaks do not
automatically indicate higher sensitivity or accuracy.

Per-library intervals and supporting QC remain separate. No consensus
intersection trims the catalogue, no IP/Input ratio becomes a significance
test, and no automatic exclusion is inferred from a generic FRiP cutoff.
Publication does not certify biological replicate independence or approve a
new paired statistical contrast.

Final native outputs, QC reports, metadata and parameters are retained
separately. Reproducible alignments and temporary data are not part of the
retained final-results panel. Gene annotation, DB, IDR, chromatin states and
atlas presentation are outside this accepted processing scope.
