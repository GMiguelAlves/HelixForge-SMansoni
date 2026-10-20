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

The former atomic support-at-least-two representation retained 367 cercarial
and 98 sporocyst segments. It has been replaced in the active results by the
nf-core/chipseq `2.1.0` merged regional catalog, minimum one calling library,
and HOMER `4.11` annotation: **1,816 cercarial and 867 sporocyst candidates**.
The ≥2-library core remains separate evidence and no longer crops the intervals
used for gene association. These are candidate regions, not independently
established biological domains. Support is between deposited libraries; their
biological independence remains unresolved. See the
[current result index](../../results/chipseq/PRJNA1090249/README.md).

**Disposition: preserve as descriptive/exploratory results with explicit
warnings; do not certify robust HP1 domain boundaries or approve DB.**
FRiP ranges from 0.201% to 3.161%, with five libraries below 0.54%, alongside
substantial differences in peak counts and the observed extension sensitivity.
FRiP alone is not used as a universal rejection cutoff, and these findings do
not establish a specific laboratory failure. They do require caution about
enrichment quality, replicate concordance and any locus-level interpretation.
The fixed extension is not retrospectively validated as a physical fragment
length by a successful run or by a nonzero consensus.

## Article-inspired follow-up comparison

A separate comparison was completed to assess the strategy used in the
publication. It does not replace the primary analysis or certify its peaks.
The comparison is motivated by the observed sensitivity above and is not a
new blinded benchmark.

The documented approach combines reads by stage for MACS2 broad calling,
whereas the primary result calls each library and then consolidates support.
The follow-up started from the preserved alignments **before the
primary MAPQ filter**, excluded alignments carrying the Bowtie2 `XS` tag,
excluded unmapped/secondary/QC-failed/supplementary records, and removed
duplicates separately within each library with `samtools rmdup -s`. There is
no extra MAPQ cutoff in this branch. An absent XS tag implements the paper's
reported uniqueness criterion; it is not an independent proof of genomic
uniqueness.

The planned target was the paper's **4.7 million reads per library**, using
Picard DownsampleSam with seed `20261018` and the `HighAccuracy` strategy.
The prespecified rule capped the common target at the minimum available count
across all eight libraries, without upsampling. A lower common depth is an
explicit deviation from the publication, not a successful reproduction of
its reported depth.

The preparation found **3,819,229 reads** after uniqueness filtering and
deduplication in the limiting library, `SRR28402818`. The prespecified rule
therefore set the common target to **3,819,229**, rather than claiming the
published 4.7 million were available in every library. Although downsampling
is probabilistic, the realized counts in this execution were exactly
**3,819,229 for each of the eight libraries**.

All three deposited IP libraries per condition are included, with the shared
Input for that condition used once. The exact subset used for the historical
pooled analysis has not been confirmed from the accessible supplementary
material. Pooling is for domain detection only; it does not create additional
biological observations or authorize differential inference.

The tools used were MACS2 `2.2.9.1` and Picard `3.4.0`; the paper does not
establish these as its original versions. Both automatic-model and separately
labeled `--nomodel --extsize 200 --shift 0` arms completed. The automatic
estimates were 72 bp for cercariae and 70 bp for sporocysts. These remain
read-length-adjacent estimates, not independently validated physical fragment
sizes. The fixed-extension arm is not a reproduction of the article's default
fragment model.

The follow-up also retains MACS2's default duplicate handling (`--keep-dup 1`),
including after pooling, and its default callpeak q-value (0.05). The latter
does not define the downstream broad score thresholds, which come from the
separate p-value score track. These are additional differences from the
primary all-tag, q-value-0.01 broad call and must not be attributed solely to
pooling or to a newer reference.

For broad calling, `bdgcmp -m ppois` produced a **negative-log10 Poisson
p-value** track from the MACS treatment pileup and control lambda, without
SPMR normalization. `bdgbroadcall` used peak-score cutoffs **2 and 5**.
These are not the primary analysis's q-value thresholds. Unreported geometry
and linking settings are explicit assumptions: link cutoff 1, minimum length
200 bp, strong-peak maximum gap 30 bp and broad linking gap 800 bp. Do not tune
these values to recover the paper's peak totals.

This is an **adaptation to SM_V10/WBPS19**, not an exact SMv9 reproduction:
reference, annotation, alignment settings, software versions and potentially
library selection/depth differ. Interval coverage, concordance and
descriptive TSS/TES patterns were examined on the current annotation. Do not compare old and
new assembly coordinates directly or claim reproduction of the original
4,898-gene metagene set. ChromstaR differential modeling and DiffBind remain
outside this follow-up's inferential scope while biological independence is
unresolved.

### Completed follow-up: results and disposition

| Condition | MACS2 model | Regions, score 2 | Regions, score 5 | Covered outer-envelope bases, score 2 |
|---|---|---:|---:|---:|
| Cercaria | automatic, d=72 | 5,080 | 0 | 6,289,399 |
| Sporocyst | automatic, d=70 | 0 | 0 | 0 |
| Cercaria | fixed 200 bp | 32,278 | 0 | 24,161,508 |
| Sporocyst | fixed 200 bp | 26,461 | 3 | 18,943,672 |

The three fixed-extension sporocyst score-5 regions cover 882 bp. **Zero called
regions is not zero local signal:** bases above score 5 were present in all
four score tracks, but most did not form intervals satisfying the declared
strong-peak length/gap requirements. This result depends on those assumptions;
the unreported original geometry was not tuned to recover the paper's totals.

For fixed-200 score-2 outer envelopes, the fraction of pooled IP reads
overlapping regions was **13.22% versus 6.91% in Input for cercariae**, and
**11.32% versus 6.19% for sporocysts**: descriptive read-fraction ratios of
**1.92 and 1.83**. These fractions use the per-library-deduplicated,
downsampled BAMs before MACS2's additional duplicate restriction after pooling.
They are not formal enrichment tests or estimates of biological uncertainty.

The fixed-200 score-2 sets overlap **33.31% of primary cercarial consensus
bases and 44.24% of primary sporocyst consensus bases**. Base Jaccard values
are only 0.000835 and 0.000505, respectively, partly reflecting the very
different sizes and granularity of the pooled outer envelopes and the small
primary support segments. Neither region counts nor Jaccard alone establish
which method is biologically correct.

Descriptive metagene profiles used RPKM coverage without fragment extension,
50 bp coverage bins, 100 bp profile bins and 5 kb flanks. The displayed
transformation is `log2(mean RPKM + 1)`, using **4,945 forward-strand
protein-coding genes** and a supplementary group of **all 9,896 protein-coding
genes** in the current annotation. The groups overlap and are not independent
samples. No zero-signal genes were removed and no biological-replicate
confidence intervals were inferred from pools. An increase near TES is
visible in both IP and Input: the aggregate shape alone does **not** establish
HP1-specific TES enrichment or reproduce the publication's conclusion.

**Processing is complete, inclusion in the atlas is accepted, and DB is not
eligible in the current analysis.** The study is not assigned an overall
exploratory-only label. Its quantitative QC and design limitations remain
attached to the results and must not be interpreted as certified enrichment
or independent biological replication. Pooling removed the computational obstacle to automatic model fitting,
but this adaptation did not establish a stable, strongly supported broad-domain
set or faithful reproduction of the published result. This does not prove
absence of biological HP1 signal, a particular laboratory failure, or a
HelixForge software defect. Primary results and this follow-up remain separate;
no parameter was selected to improve agreement with published peak counts.

## References

- [Study publication](https://doi.org/10.1590/0074-02760240075).
- [MACS3 callpeak options](https://macs3-project.github.io/MACS/docs/callpeak.html):
  model-free extension, shift, duplicate handling and broad cutoffs. The
  executed version is recorded above; the online documentation may describe
  a newer release.
