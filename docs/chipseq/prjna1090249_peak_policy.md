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
No robustness claim is made until this review is completed.

Exact boundaries and marginal loci remain sensitive to the unmeasured
fragment distribution. Findings must be labeled exploratory, particularly
where extension sensitivity or duplicate-rich enrichment affects interpretation.
Preserve the automatic-model diagnostics, commands and versions alongside the
model-free results. Do not replace the original failed-run evidence with a
successful-exit claim.

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

## References

- [Study publication](https://doi.org/10.1590/0074-02760240075).
- [MACS3 callpeak options](https://macs3-project.github.io/MACS/docs/callpeak.html):
  model-free extension, shift, duplicate handling and broad cutoffs. The
  executed version is recorded above; the online documentation may describe
  a newer release.
