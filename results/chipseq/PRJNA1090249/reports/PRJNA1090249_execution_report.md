# PRJNA1090249 — current regional results

The six original per-library MACS3 calls and matched Inputs are unchanged.
The active publication replaces atomic support-segment annotation with merged
regional candidates following nf-core/chipseq `2.1.0`, minimum one calling
library, followed by HOMER `4.11` annotation using SM_V10/WBPS19.

| Group | Regions | One-library regions | Regions called by ≥2 libraries | Associated gene IDs |
|---|---:|---:|---:|---:|
| Cercariae | 1,816 | 1,589 | 227 | 1,480 |
| Sporocysts | 867 | 793 | 74 | 771 |

The 60,657 bp cercarial and 21,643 bp sporocyst shared cores remain separate
evidence products, not the boundaries used for nearest-gene association.
Retention of the complete input union follows the minimum-one policy; it is
not a measurement of biological recall. More associated genes do not establish
greater biological accuracy or reproduction of the published study.

The source native workflow completed all 38 tasks across six groups. Direct
HOMER comparison passed for both HP1 groups; the independent H3K27me3 groups
are reported separately. Source executions predate automatic TSS ambiguity
assessment, and their explicit `not_assessed` flags preserve that limitation.
Publication only adds flattened library membership and compression; native
intervals, HOMER categories, distances and gene assignments are unchanged.

Raw-read QC, alignment, MACS3 parameters, FRiP, extension sensitivity and
article-inspired diagnostics remain available through the [result index](../README.md).
Differential Binding was not run. The reduced public package is not a new
validated native terminal run or a self-contained re-entry bundle.

Exact execution source commit, reference hashes, public/source artifact hashes,
policy and scope are recorded in each regional group manifest. The previous
atomic-consensus report and terminal validation do not certify this replacement
and have therefore been removed from the active result tree.
