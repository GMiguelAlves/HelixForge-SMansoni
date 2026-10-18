# Differential Binding review and artifact retention

Differential Binding (DB) is a planned scientific follow-up for eligible
ChIP-seq studies, not a permanently excluded analysis. Complete peak processing
and review each study promptly; eligible comparisons need not wait for other
studies. The initial consensus/IDR workflow may finish without DB while its
inputs remain protected.

## Candidates and approval

The [candidate registry](../../config/chipseq/differential_binding_candidates.tsv)
records comparisons for the five peak-producing studies. Review each mark and
biological comparison separately. The [retention policy](../../config/chipseq/differential_binding_review.json)
requires both QC and model review: mapping/complexity, Inputs, enrichment,
replicate concordance, sample identity, technical-run aggregation, compatible
references/counting units, estimable designs and normalization assumptions.
No single mapping-rate, FRiP or generic QC label automatically approves DB.

Record `ELIGIBLE`, `HOLD` or `NOT_ELIGIBLE`, the evidence, any exclusions and
their rationale. Technical libraries and shared Inputs are not additional
biological observations. Batch or pairing may be used only when supported by
metadata and an estimable model; do not infer pairs from matching replicate
numbers. PRJNA219440 and PRJNA384372 remain QC/alignment-only without a separate
design review. The PRJNA816041 technical benchmark is not an automatic DB study.

## HP1: PRJNA1090249

The [prospective plan](../../config/chipseq/PRJNA1090249/db_plan.json) registers
**cercaria versus sporocyst**, three biological HP1 IP replicates per stage,
with the initial unpaired model `~ condition`. Positive log2 fold change means
higher relative HP1 enrichment in cercaria. The two shared Inputs remain
controls, not columns in the inferential model. Their lack of independent
within-stage replication is a limitation to record, not a reason to duplicate
them computationally.

The candidate peak universe is the union of accepted condition-specific
consensus sets. Count the same universe in all six IP BAMs using single-end
read units; do not infer DB simply from presence/absence in called peak sets.
QC failure or loss of a planned replicate places the decision on hold rather
than silently reducing the design. A compatible executable statistical
specification, count filtering, normalization review and R/counting runtime
preflight are required before execution. The plan is **not** an executable
`chipseq_db_spec` and does not enable DB in a running workflow.

Standard count normalization requires review when widespread changes are
plausible; relative DB must not be presented as calibrated absolute occupancy.
See the [DESeq2 guidance on size-factor assumptions](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#control-features-for-estimating-size-factors)
and [non-estimable designs](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#model-matrix-not-full-rank).

## Retention before cleanup

Preserve final **BAM/BAI per biological replicate**, control BAM/BAI, producer
manifests/checksums, normalized metadata, control maps, the peak plan, reference
identity, peaks/consensus manifests, QC, parameters and versions until DB is
resolved. Do not replace replicate BAMs with pooled BAMs. Preserve existing
workdir dependencies until an independent retained-artifact input inventory is
verified. Missing review or `HOLD` does not permit BAM cleanup.

After an explicit decision not to run DB, or successful DB and verified
persistence of its results/audit, release disposable intermediates under the
normal storage policy. Finish the review/DB cycle study by study; if an
unresolved hold would exceed storage capacity, request a retention decision
rather than deleting the evidence or staging the next heavy study blindly.

The current HelixForge DB subworkflow accepts declared BAM, consensus and
metadata channels. The `post_consensus` entry point performs annotation/tracks/
reporting, **not DB**. Re-entry for DB must use the existing DB subworkflow
with verified retained artifacts, or a verified cache-preserving execution;
do not assume that passing a different run mode guarantees upstream cache reuse.
