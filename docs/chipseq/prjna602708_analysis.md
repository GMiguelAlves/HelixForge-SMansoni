# PRJNA602708: ChIP-seq analysis and mapping sensitivity

## Study and reference

This study contributes paired-end H3K4me3 and H3K27me3 libraries from
schistosomula exposed to plasma collected from rhesus macaques or infected
golden hamsters. The five technical condition identifiers are `hm7w`, `rh0w`,
`rh8w`, `rh10w`, and `rh12w`, with two deposited IP replicates per mark and
condition. The original descriptions are retained in
[`design.tsv`](../../metadata/chipseq/PRJNA602708/design.tsv).

The primary publication is [Amaral et al., 2021](https://doi.org/10.1038/s41467-021-26497-0).
All reanalyses use `SM_V10/WBPS19`, reference identifier
`Schistosoma_mansoni_SM_V10_WBPS19`. They are not an exact reconstruction of
the publication's reference, historical AQUAS configuration or peak catalogue.

## Two separate result layers

The study retains both mapping policies, rather than silently selecting one
because it returns more peaks or greater replicate overlap:

| Layer | Alignment selection | Interpretation |
| --- | --- | --- |
| `mapq30` | MAPQ >=30; exclude flags2308 | More restrictive mapping-confidence selection, not a complete or independently validated occupancy catalogue |
| `mapq0` | No minimum MAPQ; exclude flags2308 | Sensitivity layer including reported low-confidence primary mappings, not proof of unique localization |

Both exclude unmapped, secondary and supplementary records. The second layer
does not recover alternative alignments that Bowtie2 did not report.
Peak IDs and gene associations must always be qualified by their layer;
combining the two layers would mix different input and background populations.

The remaining study-specific choices are identical: samtools1.23.1 duplicate
removal within each IP/Input library, mitochondrial read/mate exclusion before
peak calling, MACS3 3.0.4 `BAMPE`, `p=0.01`, `--keep-dup all` after upstream
deduplication, and effective genome size391422299. That size is an assembly-span
proxy, not a validated mappability estimate. The broad weak-region threshold
remains `p=0.1`; neither threshold is an empirical FDR guarantee.

## Controls and experimental design

Two Inputs are shared: replicate1 IPs use `SRR10961878`; replicate2 IPs use
`SRR10961877`, across conditions and marks. This association is the project's
registered reconstruction, not a recovered command-level IP/Input map from
the authors. The two no-antibody/beads libraries remain QC negative controls
and are not substituted for Inputs.

Matching replicate numbers identify this control assignment; they do not
establish statistical pairing of biological observations across conditions.
No cross-condition matched biological unit is established by this metadata.
Potential treatment/time comparisons would therefore require an independently
supported design, not a paired model inferred from library names.

Differential Binding remains **suspended**. Technical completion, two deposited
replicate labels or a larger peak catalogue do not establish eligibility for
inferential DB. Peak presence/absence must not be described as a significant
between-condition difference.

## Measured sensitivity

The table reports sums of individual library calls, not unique regions,
consensus sizes or associated genes.

| Policy comparison | H3K4me3 calls | H3K27me3 calls |
| --- | ---: | ---: |
| MAPQ30, individual Inputs | 22,672 | 377,974 |
| MAPQ20, individual Inputs | 23,853 | 421,164 |
| No minimum MAPQ, individual Inputs | 14,532 | 429,300 |
| MAPQ30, pooled Inputs, separate sensitivity only | 99,358 | 827,933 |

H3K4me3 replicate coverage Jaccard changes from 0.07–0.45% at MAPQ30 to
6.36–31.98% without a minimum MAPQ. However, only 7–35% of the MAPQ30 group-union
coverage is retained in the no-minimum layer. The effect is spatial, not just
a change in the number of peaks. Higher concordance does not distinguish
restored signal in repetitive sequence from shared ambiguous signal or artifacts.

H3K27me3 responds heterogeneously: no-minimum Jaccard improves in four groups
and decreases in `rh12w`. The pooled-Input sensitivity changes control depth,
MACS scaling and local background; it is not adopted or mixed into either layer.

Direct MACS3 execution reproduced the original twenty HF calls byte-for-byte
across 100 native scientific output files. This establishes wrapper equivalence
for those inputs/arguments, not biological correctness or equivalence to the
published catalogue. The controlled investigation is complete with limitations;
no configuration is labelled a proven reproduction of the article.

## Regional representation and gene associations

The regional policy is the pinned nf-core/chipseq 2.1.0 implementation with
minimum-one-library support. Complete merged candidate regions and their
per-library membership/coverage are retained. The shared core is a separate
spatial evidence product, not a rule that clips the annotation catalogue.

HOMER 4.11 with the frozen GTF supplies gene assignments. Each mapping layer
has its own annotation and individual CPM tracks; historical annotations were
not reused after peak coordinates changed. Candidate regions, library support and
assignment-ambiguity fields are evidence, not proof of occupancy or regulation.

The [public results](../../results/chipseq/PRJNA602708/README.md) include twenty
regional catalogs, forty individual IP peaksets, forty-eight policy/mark QC
records, four native HTML reports and a track inventory. All 464 native
postprocessing tasks completed. Shared Inputs appear in both mark-specific
QC sets; these are not forty-eight distinct biological observations.
Heavy BAMs/BAIs and BigWigs are retained outside Git. The public copy is not
a self-contained re-entry package.

IDR from the original sparse q-value calls was not evaluable. It must not be
reported as successful, as zero reproducible regions, or as an IDR result for
the newer p-value layers. No new IDR result is claimed for this delivery.
