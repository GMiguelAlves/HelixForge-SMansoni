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

### Required final design review for every study

Every final study review must document:

- **Pairing decision and rationale:** paired or unpaired, with the publication
  URL/DOI, methods section and supporting design evidence. Pairing requires a
  documented matched biological unit or blocking relationship; replicate
  numbers alone do not establish one.
- **Experimental unit and independence:** map each sequenced library to its
  source individual, independently collected pool or preparation. Identify
  aliquots, repeated sampling and technical libraries explicitly. Distinct
  sample IDs or concordant profiles do not prove biological independence.
- **Confounding and interpretation:** identify factors that cannot be separated
  and label the contrast accordingly. A full-rank group-only model does not
  make a confounded biological effect identifiable.
- **Eligibility:** record QC evidence, the proposed model and the decision.
  Unresolved biological independence means `HOLD` for inferential DB, even
  when technical QC passes. Descriptive peak processing can remain available
  with its limitations; it is not certification of biological replication.

## HP1: PRJNA1090249

The [HP1 peak-calling policy](prjna1090249_peak_policy.md) documents the
model-free fragment-extension assumption and its planned sensitivity review.
This descriptive processing does not lift the inferential DB hold below.

The [prospective plan](../../config/chipseq/PRJNA1090249/db_plan.json) registers
**cercariae SmVEN versus sporocysts SmDFO**, with three deposited HP1 IP
libraries per group. The [publication's methods, page 2](https://pdfs.semanticscholar.org/91d7/7351164c927b29796157e2408e27c8bffa9e.pdf)
describe sporocyst dissection from snails carrying SmDFO and light-induced
cercarial shedding using SmVEN. They do not establish matched biological units
across stages. The proposed comparison is therefore **unpaired**, using
`~ condition` only if biological independence and the remaining gates are
confirmed. Existing technical condition IDs remain unchanged.

**Stage and parasite strain are completely confounded.** Any DB result would
describe differences between these combined groups, not an isolated
developmental effect. Adding strain and stage to the design or applying batch
correction cannot separate their contributions. Positive log2 fold change
would mean higher relative HP1 enrichment in cercariae SmVEN than sporocysts
SmDFO. The two shared Inputs remain
controls, not columns in the inferential model. Their lack of independent
within-stage replication is a limitation to record, not a reason to duplicate
them computationally.

**Biological independence is not yet confirmed.** The methods describe two
sporocyst samples per snail without linking these preparations to the deposited
libraries. Cercarial aliquots likewise require source tracing. The three
replicate labels in each group must not be presented as three established
independent biological observations. DB is on `HOLD` until library-to-source
mapping resolves this uncertainty. If libraries share a biological source,
review aggregation and the resulting number of independent units before any
inference; do not silently retain a nominal three-versus-three model.

The candidate peak universe is the union of accepted condition-specific
consensus sets. If approved, count the same universe in the retained IP BAMs using single-end
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

Preserve final **BAM/BAI per library and confirmed biological unit**, control BAM/BAI, producer
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
