# MACS3 fragment-model sensitivity — PRJNA1090249

No production parameter, original BAM or peak set was changed. This test only
assessed whether limiting duplicate tag coordinates, at the original modeling
settings, restored an estimable MACS3 model. It did not.

| Library | Paired regions, all tags | Selected d | Paired regions, one tag per position/strand | Selected d |
|---|---:|---:|---:|---:|
| Cercaria HP1 r1 | 142 | 72 bp | 88 | not estimable |
| Cercaria HP1 r2 | 119 | 71 bp | 80 | not estimable |
| Cercaria HP1 r3 | 118 | 71 bp | 77 | not estimable |
| Sporocyst HP1 r1 | 81 | not estimable | 45 | not estimable |
| Sporocyst HP1 r2 | 74 | not estimable | 41 | not estimable |
| Sporocyst HP1 r3 | 79 | not estimable | 57 | not estimable |

MACS3 requires at least 100 paired regions for this model. All six technical
diagnostics completed, but completion is not a successful fragment estimate.
The model outputs and logs, not only process exit codes, determine estimability.

MACS3 3.0.4 `predictd` used genome size 391422299, mfold 5–50, bandwidth 300 and
minimum d 20. The sensitivity arm used `macs3 filterdup --keep-dup 1` to generate
temporary BED representations. These were removed after diagnostic completion;
original BAMs remained unchanged. Coordinate-based limiting is not identical
to historical `samtools markdup` and does not reproduce the full legacy run.

No mfold sweep, new peak calling, nomodel or fixed extsize was applied. These
tests do not certify the retained 71–72-bp estimates, prove the libraries are
unusable, or identify the cause of the difference from a successful older run.
Final fragment length remains **INCONCLUSIVE** pending further evidence or an
explicit justified modeling decision.
