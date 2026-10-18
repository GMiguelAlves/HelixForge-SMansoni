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

## References

- [Study publication](https://doi.org/10.1590/0074-02760240075).
- [MACS3 callpeak options](https://macs3-project.github.io/MACS/docs/callpeak.html):
  model-free extension, shift, duplicate handling and broad cutoffs. The
  executed version is recorded above; the online documentation may describe
  a newer release.
