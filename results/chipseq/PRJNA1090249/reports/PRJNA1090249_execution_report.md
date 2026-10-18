# PRJNA1090249 — final execution review

## Scope and completed processing

Eight single-end libraries passed the download integrity gate: six HP1 IPs
and two stage-specific Inputs. The native workflow completed QC, trimming,
Bowtie2 alignment, BAM processing, MACS3 broad calls, FRiP, consensus,
annotation, coverage tracks, reporting and terminal-manifest validation.
The terminal validation passed both schema and semantic checks. The final
model-free continuation reused 97 upstream tasks rather than realigning.

HelixForge `v1.2.0`, Nextflow `25.10.7`, Java `21` and the certified Conda
runtimes were used. SM_V10/WBPS19 reference identity and exact provenance are
available in the public reference and execution manifests.

## Primary peak calling

MACS3 `3.0.4`, broad mode, q-value `0.01`, broad cutoff `0.1`,
`--nomodel --extsize 200 --shift 0`, `--keep-dup all`; MAPQ >= 30,
excluded SAM flags 2308. The effective-size value 391,422,299 is an assembly
span proxy, not a measured mappable-genome size. No validated blacklist was
available. The fixed extension is not an experimentally measured fragment size.

| IP library | Peaks | FRiP (%) |
|---|---:|---:|
| cercaria_hp1_r1 | 618 | 0.3926 |
| cercaria_hp1_r2 | 267 | 0.2236 |
| cercaria_hp1_r3 | 1,240 | 0.5364 |
| sporocyst_hp1_r1 | 500 | 0.2422 |
| sporocyst_hp1_r2 | 407 | 0.2010 |
| sporocyst_hp1_r3 | 48 | 3.1607 |

Consensus with support >= 2 retained 367 cercarial and 98 sporocyst atomic
support segments. Their covered lengths are 60,657 and 21,643 bp,
respectively. Segment counts must not be interpreted as counts of complete
biological broad domains.

## Follow-up and interpretation

Extension sensitivity at 150 and 250 bp produced base Jaccard values
0.33558–0.69351 relative to 200 bp. No setting was chosen to improve agreement
with published counts. The article-inspired pooled MACS2 comparison used a
common realized depth of 3,819,229 reads per library, not the publication's
4.7 million. The fixed-200 arm produced 32,278 and 26,461 score-2 regions for
cercariae and sporocysts, but zero and three score-5 regions. Unreported
original interval geometry, reference and software differences prevent a
claim of exact reproduction. Neither divergence nor zero qualifying domains
proves absence of local HP1 signal or a particular laboratory failure.

IP/Input read-fraction ratios in fixed-200 score-2 outer envelopes were 1.92
and 1.83. The TSS/TES profiles show TES-associated coverage in Input as well
as IP. These are descriptive measurements, not formal enrichment tests or
biological-replicate confidence estimates. The complete methods, figures
and tables are linked from the [result index](../README.md).

## Final design and use decision

Processing is closed and atlas inclusion is accepted. No global
exploratory-only classification is assigned to the study. The quantitative
QC findings remain available and must accompany interpretation.

The potential group comparison would be **unpaired cercariae SmVEN versus
sporocysts SmDFO**, not an isolated developmental-stage effect. The source
methods do not establish matched units across stages; stage and strain are
fully confounded. Within-group biological independence is unresolved. Inputs
are controls, not additional biological replicates. Therefore DB is
`NOT_ELIGIBLE_IN_CURRENT_ANALYSIS`, and no inferential DB results are claimed.
Reassessment requires new source or experimental evidence rather than more
parameter tuning.

## Retention

Scientific results, figures, reports, annotation, terminal provenance and
coverage tracks are retained. Final BAM/BAI files for all eight libraries
are preserved privately with verified checksums. Private audit is archived
separately; reproducible intermediates are removed only after preservation
is verified. No private filesystem paths are part of this public review.

[Study publication and source methods](https://doi.org/10.1590/0074-02760240075).
