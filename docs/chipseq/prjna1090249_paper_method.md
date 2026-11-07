# PRJNA1090249: article-guided HP1 reanalysis specification

## Source and scope

This prospective processing specification follows **da Trindade et al. (2025)**,
*Heterochromatin protein 1 (HP1) of Schistosoma mansoni: non-canonical chromatin
landscape and oviposition effects*, DOI
[10.1590/0074-02760240075](https://doi.org/10.1590/0074-02760240075), especially
the **Bioinformatic analyses** methods and the results on broad HP1 regions.
The [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11961034/) and its
Supplementary Table I are the methodological sources.

The scope is six single-end HP1 IP libraries and two stage-specific Inputs:

| Group | IP libraries | Input |
|---|---|---|
| Cercariae | SRR28402826, SRR28402825, SRR28402822 | SRR28402821 |
| Sporocysts | SRR28402820, SRR28402819, SRR28402818 | SRR28402817 |

The [curated design](../../metadata/chipseq/PRJNA1090249/design.tsv) preserves
library identity. Drosophila antibody controls are not downloaded as part of
this S. mansoni reanalysis. Existing accepted HF results remain separate until
the new outputs have been reviewed; this document does not declare a new result.

## Published processing

The article describes 75 bp single-end NextSeq 550 sequencing, FastQC/MultiQC,
Trim Galore adapter removal, sensitive end-to-end Bowtie2 alignment to **SMv9 /
WBPS17**, retention of unique alignments using the `XS` tag, SAMtools `rmdup`,
and Picard downsampling to **4.7 million aligned unique reads per library**.

For broad regions, the authors combine uniquely aligned ChIP reads **by stage**
and apply MACS2 followed by `bdgbroadcall`, with score cutoffs **2 and 5**.
The article interprets these as p-value scores 10^-2 and 10^-5, not MACS q-values.
Its separate differential analyses use ChromstaR (1,000 bp bins, 500 bp steps)
and DiffBind; those inferential analyses are not prerequisites for producing
the stage-level region catalogues requested here.

Published region counts are 7,485/222 for cercariae and 6,406/89 for sporocysts
at the two thresholds. They are comparison references, **not optimization
targets or acceptance gates**.

## Explicit implementation adaptations

The reanalysis keeps **SM_V10/WBPS19** and existing installed tool versions.
Their exact runtime identities must be recorded; no historical-version matching
or replacement by MACS3 is implied. Reuse the compatible Bowtie2 index.

Available version probes identify FastQC 0.12.1, Bowtie2 2.5.5, SAMtools 1.23.1,
MACS2 2.2.9.1 and Picard 3.4.0. `samtools rmdup -s` remains available in this
installation. These are reconstruction versions, not versions established for
the historical Galaxy execution.

1. **QC/trimming:** use the registered
   [ChIP-specific policy](reference_and_trimming.md): SE Trim Galore, Q20,
   minimum length 25 bp and adapter autodetection, with raw/post-trim QC.
   These numerical choices are explicit reconstruction settings because the
   article does not specify them; they are not copied from RNA-seq. Inspect
   tagmentation adapters and NextSeq artifacts rather than silently adding flags.
2. **Alignment/filtering:** sensitive end-to-end Bowtie2; discard unmapped,
   secondary, supplementary and QC-failed records and alignments carrying `XS`.
   Do not add the old HF MAPQ >=30 rule to this article-guided branch. The flag
   exclusions operationalize primary usable alignments; the paper does not
   provide their command. Absence of `XS` is not independent proof of unique
   genomic localization.
3. **Duplicates/depth:** use single-end `samtools rmdup -s` and count available
   reads per library before Picard downsampling. The target remains 4,700,000;
   use Picard's HighAccuracy strategy with reconstruction seed 20261107 and
   record the actual retained counts. These sampling details are not reported
   by the article. If a library
   has fewer reads, report the shortfall before choosing a depth exception.
   Do not upsample or silently reduce all libraries to the smallest count.
4. **Pooling:** retain per-library QC and include all three IP libraries per
   stage in the pooled call. Use that stage's shared Input once. Pooling is
   for region detection, not extra biological replication. The paper does not
   provide an exact executable map of every pooled file.
5. **MACS2:** retain automatic fragment modeling initially, explicit assembly-span
   proxy `g=391422299` for SM_V10, and record the installed defaults, including
   internal duplicate handling. This proxy is not a measured effective/mappable
   genome size. Do not automatically substitute the earlier fixed 200 bp model
   or `--keep-dup all`. Model failure requires review, not an invisible fallback.
6. **Broad score:** generate a documented negative-log10 Poisson p-score track
   with MACS2 `bdgcmp -m ppois`, then call `bdgbroadcall` at 2 and 5. The intermediate
   score-generation command is an explicit reconstruction step, not a recovered
   historical command. Freeze the installed tool's unreported length/gap/link
   defaults before calls; do not reuse earlier sensitivity overrides as if the
   article specified them. Do not relabel a treatment pileup, fold-enrichment
   track or q-score as this p-score.

All native regions and both thresholds remain separately identified. No HF
shared-core intersection or support-at-least-two rule recuts the pooled regions.
Optional HOMER annotation on WBPS19 is an additional project output, not the
paper's original gene-association procedure. It must preserve threshold,
stage, library provenance, native coordinates and association uncertainty.
Pool membership is not evidence that each region is supported by every
constituent library; do not label it as per-region replicate support.

## Interpretation and earlier evidence

The [earlier HP1 policy and comparison](prjna1090249_peak_policy.md) already
record an article-inspired pooled sensitivity analysis. It used a lower common
depth and explicit fragment/geometry assumptions; it is not the new protocol
and does not establish reproduction. Its numerical limitations should inform
review, not silently dictate new parameter choices.

The article uses cercariae from SmVEN and sporocysts from SmDFO. Stage and strain
are therefore confounded. Replicate labels do not establish cross-stage pairing
or prove independent biological sources. This processing does not authorize DB,
an isolated developmental-effect claim or an exact historical replication.

Retain final region/annotation tables, figures, QC metrics and compact provenance.
Do **not** retain BAM/BAI, FASTQs, workdirs or heavy fitted objects as final
archives. Future alignment-dependent analyses would require reprocessing.
