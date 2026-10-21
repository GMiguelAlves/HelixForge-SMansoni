# PRJNA312093 — Differential Binding eligibility

Peak processing is complete for H3K27me3 and H3K4me3. The regional catalogs,
library support, HOMER associations, H3K4me3 IDR results and coverage tracks
remain usable as measured evidence. **Inferential Differential Binding (DB)
is on HOLD**, independently of atlas inclusion. See the
[machine-readable review](../../config/chipseq/PRJNA312093/db_review.json).

## Experimental units and pairing

The proposed contrasts are male versus female, separately for each stage and
mark. They are **unpaired**: males and females are distinct parasite pools;
matching replicate labels do not establish matched samples. Inputs are
controls, never additional observations in the model.

[Picard et al. 2016](https://doi.org/10.1371/journal.pntd.0004930), Methods,
describes unisexual infections, pooled parasites and biological replication;
its ChIP description reports two male and three female replicates. The
[2019 H3K4me3 methods](https://doi.org/10.1093/gbe/evz133) explicitly describe
two male biological replicates and distinguish sexually immature unpaired
worms from mature paired adults. The selected ENA records identify two female
libraries per stage, but do not resolve their exact independent source pools
against those published designs. All selected runs refer to a generic
[BioSample](https://www.ebi.ac.uk/ena/browser/view/SAMN04493444). That shared
identifier does **not** prove that the libraries are technical replicates;
neither does it verify biological independence.

## Decisions by comparison

| Mark and stage | Decision | Reason |
|---|---|---|
| H3K4me3, adult | HOLD | Male IPs are HiSeq 2000; female IPs are HiSeq 2500. Sex and instrument cannot be separated. Source-pool mapping also remains unresolved. |
| H3K27me3, adult | HOLD | The same complete sex–instrument confounding and source-pool issue apply. |
| H3K4me3, cercaria | HOLD | All four IPs use HiSeq 2500 and enrichment metrics are favorable; exact source-pool independence must still be established. |
| H3K27me3, cercaria | HOLD | Source-pool mapping is unresolved. One male IP is HiSeq 2000; the other three IPs are HiSeq 2500. An adjusted four-sample model has only one residual degree of freedom. |

Instrument assignments come from the
[frozen ENA records](../../provenance/source_package/frozen/chipseq/PRJNA312093/ena_read_run.tsv).
For each adult comparison, `~ instrument + sex` has rank 2 for 3 columns:
adding a batch term cannot repair complete confounding. A sex-only fit could
describe a combined group difference, but would not identify an isolated sex
effect. No such alternative inference is silently substituted.

## QC and requirements for reconsideration

H3K4me3 selected BAMs contain 9.8–26.7 million reads per IP. Per-library FRiP
is 0.278–0.571. IDR at 0.05 retains 4,348 female-adult, 4,461
female-cercaria, 4,719 male-adult and 4,651 male-cercaria regions. These are
positive enrichment/reproducibility observations, not proof of independent
biological sampling or of every minimum-one regional candidate.

H3K27me3 FRiP ranges from 0.021 to 0.213. Adult female libraries differ
substantially (0.194 versus 0.021), as do adult males (0.092 versus 0.030).
This is an additional quality-review concern, not a newly invented universal
FRiP failure threshold for broad marks. Retaining a regional union does not
establish concordant quantitative behavior between those libraries.

Duplicates were retained. A zero duplicate flag in an unmarked BAM must not
be interpreted as high library complexity. If DB becomes eligible, library
complexity, Input behavior, a common counting universe, sample concordance,
count filtering and normalization assumptions require a complete review.
Without calibration, DB would describe relative enrichment, not absolute
chromatin occupancy.

Reconsideration requires a library-to-independent-source map and a defensible
handling of the instrument-confounded comparisons, followed by model and
runtime validation. Final per-library IP/Input BAMs, indexes, metadata and
peak evidence remain protected. Disposable workdirs may be removed only after
their needed artifacts have verified persistent copies and a validated
retained-artifact input inventory. No DB results or statistical sex effects
are claimed by the current peak publication.
