# HP1 peak-calling policy: PRJNA1090249

## Scope and rationale

HP1 domains are processed descriptively against SM_V10/WBPS19 with matched
Inputs. A successful peak-calling command does not by itself establish
library quality, biological replication or eligibility for Differential
Binding (DB). The [DB review](differential_binding_review.md) remains separate.

The initial MACS3 single-end model could not be built for any of the three
sporocyst IP libraries: 74–81 paired strand-enrichment regions were detected,
below the 100 required by the model. Cercarial libraries produced estimates
of 71–72 bp, close to the observed modal read length of 76 bp. Exploratory
strand-correlation diagnostics did not establish a reliable, shared physical
fragment length. These observations are not proof of an irrecoverably poor
library or of a specific library-preparation failure.

A model-only sensitivity analysis limited tags to one per genomic
coordinate and strand using MACS3 `filterdup`. None of the six IP libraries
then provided enough paired regions for the default model. This does not
reproduce every aspect of `samtools markdup`, but it does not support
duplicate restriction alone as a solution to the modeling failure.

No sample-specific BioAnalyzer profiles were recovered from the examined
publication, supplement or public sequencing metadata. Their absence from
those sources does not imply they were never generated.

## Explicit model-free assumption

Before examining new model-free peak results, the primary configuration was
fixed as:

```text
MACS3 3.0.4
format: BAM (single-end)
peak type: broad
additional arguments: --nomodel --extsize 200 --shift 0
q-value: 0.01
broad cutoff: 0.1
genome-size argument: 391422299 (assembly-span proxy)
duplicate policy: all
```

**200 bp is an operational assumption, not an experimentally measured or
independently validated fragment length.** It follows the conventional
fixed-extension setting supported by MACS, not an optimization against the
number of peaks or their significance. MACS3 extends single-end tags in the
5′→3′ direction under `--nomodel`; shift zero retains their original starts.
This is not the centered-cut-site convention used for some accessibility assays.

The extension assumption is applied consistently to all six IPs, including
cercariae; successful but suspect automatic estimates are not mixed with
fixed-extension sporocyst results. Matched Inputs, mapping filters, reference,
significance thresholds and broad-peak calling remain unchanged. Existing
BAMs are not deduplicated or overwritten. Retaining the existing duplicate
policy isolates the model-free change; it is not evidence that every
duplicate alignment is biological signal.

## Bounded sensitivity review

The prespecified alternative extensions are **150 and 250 bp**, applied to
all six IP libraries with otherwise identical settings. Report peak counts,
covered bases and base/interval concordance with the 200 bp primary results.
These are sensitivity descriptions, not a competition to select the setting
producing more peaks, more significance or a preferred biological conclusion.
The completed review below does **not** establish extension-invariant peaks.

Exact boundaries and marginal loci remain sensitive to the unmeasured
fragment distribution. Findings must be labeled exploratory, particularly
where extension sensitivity or duplicate-rich enrichment affects interpretation.
Preserve the automatic-model diagnostics, commands and versions alongside the
model-free results. Do not replace the original failed-run evidence with a
successful-exit claim.

## Completed descriptive processing and sensitivity results

All six primary model-free calls completed, followed by replicate-support
consensus, annotation, tracks, an HTML report and a validated terminal manifest.
All twelve alternative calls and their comparison also completed using the
same native MACS3 provider. Treatment and matched Input BAM checksums were
verified against the primary execution. The primary setting remains 200 bp;
no setting was selected from the resulting peak counts or significance.

| Library | Peaks, 150 bp | Peaks, 200 bp | Peaks, 250 bp | FRiP, 200 bp | Base Jaccard, 150 vs 200 | Base Jaccard, 250 vs 200 |
|---|---:|---:|---:|---:|---:|---:|
| Cercaria r1 (`SRR28402826`) | 659 | 618 | 741 | 0.393% | 0.393 | 0.442 |
| Cercaria r2 (`SRR28402825`) | 255 | 267 | 262 | 0.224% | 0.389 | 0.442 |
| Cercaria r3 (`SRR28402822`) | 1521 | 1240 | 967 | 0.536% | 0.392 | 0.468 |
| Sporocyst r1 (`SRR28402820`) | 738 | 500 | 379 | 0.242% | 0.363 | 0.461 |
| Sporocyst r2 (`SRR28402819`) | 629 | 407 | 362 | 0.201% | 0.336 | 0.390 |
| Sporocyst r3 (`SRR28402818`) | 40 | 48 | 26 | 3.161% | 0.682 | 0.694 |

Base Jaccard is the intersection divided by the union of covered bases on
matching chromosomes, using zero-based, half-open coordinates. Overlapping
intervals within each set are merged before counting bases. It is not the
fraction of identical peak records. Across these twelve comparisons, the
alternatives retain **49.3–74.8%** of the bases called at 200 bp. Thus the
extension assumption materially affects the called sets, not only their
counts. No retrospective numerical PASS threshold is imposed on these metrics.

At 200 bp, the support-at-least-two consensus retains **367 intervals for
cercariae and 98 for sporocysts**. These are the consensus interval outputs,
not a count of independently established biological domains. Support is
between deposited libraries; their biological independence remains unresolved.

**Disposition: preserve as descriptive/exploratory results with explicit
warnings; do not certify robust HP1 domain boundaries or approve DB.**
FRiP ranges from 0.201% to 3.161%, with five libraries below 0.54%, alongside
substantial differences in peak counts and the observed extension sensitivity.
FRiP alone is not used as a universal rejection cutoff, and these findings do
not establish a specific laboratory failure. They do require caution about
enrichment quality, replicate concordance and any locus-level interpretation.
The fixed extension is not retrospectively validated as a physical fragment
length by a successful run or by a nonzero consensus.

## Historical analysis and inferential limits

The [legacy repository](https://github.com/GMiguelAlves/chipseq_pipeline),
examined at commit `7eac25b1294335d4c63d94cd4c497076890ed515`, supports
duplicate removal and extra MACS options. Its defaults do not establish the
actual settings of a previous user execution. This analysis does not claim
numerical equivalence to that execution or to the publication's earlier
reference and processing strategy.

DB remains on hold pending confirmation of library-to-biological-source
mapping. Cercariae SmVEN versus sporocysts SmDFO is an **unpaired combined-group
comparison**: stage and strain are completely confounded. Model-free peak
calling neither resolves that confounding nor proves three independent
biological replicates per group. Replicate-support consensus is descriptive
library support until independence is established.

## Article-inspired follow-up comparison

A separate comparison is being prepared to assess the strategy used in the
publication. It does not replace the primary analysis or certify its peaks.
The comparison is motivated by the observed sensitivity above and is not a
new blinded benchmark.

The documented approach combines reads by stage for MACS2 broad calling,
whereas the primary result calls each library and then consolidates support.
The follow-up therefore starts from the preserved alignments **before the
primary MAPQ filter**, excludes alignments carrying the Bowtie2 `XS` tag,
excludes unmapped/secondary/QC-failed/supplementary records, and removes
duplicates separately within each library with `samtools rmdup -s`. There is
no extra MAPQ cutoff in this branch. An absent XS tag implements the paper's
reported uniqueness criterion; it is not an independent proof of genomic
uniqueness.

The target is the paper's **4.7 million reads per library**, using Picard
DownsampleSam with seed `20261018` and the `HighAccuracy` strategy. Realized
counts must be recorded. If any of the eight libraries has insufficient
post-deduplication reads, the common target becomes the minimum available
count, capped at 4.7 million. Never duplicate reads to fill a library. A lower
common depth is an explicit deviation from the publication, not a successful
reproduction of its reported depth.

All three deposited IP libraries per condition are included, with the shared
Input for that condition used once. The exact subset used for the historical
pooled analysis has not been confirmed from the accessible supplementary
material. Pooling is for domain detection only; it does not create additional
biological observations or authorize differential inference.

The intended pinned tools are MACS2 `2.2.9.1` and Picard `3.4.0`; the paper
does not establish these as its original versions. MACS2 default fragment
modeling will be attempted and its diagnostics retained. If modeling is not
usable, a separately labeled `--nomodel --extsize 200 --shift 0` arm may be
evaluated; it must not be reported as successful reproduction of the article's
default model.

For broad calling, `bdgcmp -m ppois` will produce a **negative-log10 Poisson
p-value** track from the MACS treatment pileup and control lambda, without
SPMR normalization. `bdgbroadcall` will use peak-score cutoffs **2 and 5**.
These are not the primary analysis's q-value thresholds. Unreported geometry
and linking settings are explicit assumptions: link cutoff 1, minimum length
200 bp, strong-peak maximum gap 30 bp and broad linking gap 800 bp. Do not tune
these values to recover the paper's peak totals.

This is an **adaptation to SM_V10/WBPS19**, not an exact SMv9 reproduction:
reference, annotation, alignment settings, software versions and potentially
library selection/depth differ. Compare interval coverage, concordance and
descriptive TSS/TES patterns on the current annotation. Do not compare old and
new assembly coordinates directly or claim reproduction of the original
4,898-gene metagene set. ChromstaR differential modeling and DiffBind remain
outside this follow-up's inferential scope while biological independence is
unresolved. Results of this follow-up are not yet available.

## References

- [Study publication](https://doi.org/10.1590/0074-02760240075).
- [MACS3 callpeak options](https://macs3-project.github.io/MACS/docs/callpeak.html):
  model-free extension, shift, duplicate handling and broad cutoffs. The
  executed version is recorded above; the online documentation may describe
  a newer release.
