# PRJNA312093 — H3K27me3 regional publication

Four groups completed regional consolidation and HOMER annotation using the
preserved per-library MACS3 peaksets. No download, alignment or MACS3 call was
repeated for this replacement. The active policy is nf-core/chipseq `2.1.0`,
minimum one calling library, with the ≥2-library shared core retained separately.

| Group | Regions | One-library regions | Regions called by both libraries | Associated gene IDs |
|---|---:|---:|---:|---:|
| Female adult | 75,155 | 71,400 | 3,755 | 6,120 |
| Female cercaria | 29,089 | 20,537 | 8,552 | 4,326 |
| Male adult | 39,845 | 34,165 | 5,680 | 5,031 |
| Male cercaria | 44,179 | 29,716 | 14,463 | 5,455 |

The complete regional intervals, per-library covered bases/fractions, calling
library IDs, support categories, original HOMER assignments and shared cores
are published. Publication adds compression and flattened membership fields;
it does not change peak boundaries, HOMER assignments, categories or distances.

The source workflow finished all 38 tasks across six groups, including the
two separately published HP1 groups. The strict direct-tool comparison passed
for the female H3K27me3 groups but detected one assignment difference in each
male group. The affected genes share the same TSS. That failure is retained in
the publication record, not reclassified as a strict pass. Automatic ambiguity
assessment was not available in this source execution and remains explicitly
unassessed. More regions or associations do not prove improved biological
accuracy or equivalence to a publication's processing strategy.

## H3K4me3 completion

The narrow arm completed 226 native tasks: 82 QC/alignment/BAM, 51
MACS3/Peak QC/regional consolidation, nine IDR and 84 post-consensus tasks.
The existing Bowtie2 index and matched Input BAMs were reused. Per-library
FRiP is 0.278–0.571, and four groups retain 4,348–4,719 IDR regions.
Native HOMER annotation, 20 coverage tracks, the continuation HTML report and
the terminal manifest were produced. Schema, semantic and tracked-file
validation passed. These are processing/contract checks, not a new independent
biological benchmark or a direct-tool comparison.

The H3K4me3 source commit and assessed/unassessed ambiguity states are recorded
in its [execution summary](../narrow/manifests/native_execution_summary.json).
The earlier H3K27me3 comparison outcome above remains unchanged. There is no
new combined native terminal report for both marks.

DB is [on HOLD](../../../../docs/chipseq/prjna312093_differential_binding.md):
adult sex–instrument confounding, unresolved selected female source-pool
mapping and comparison-specific QC/model concerns prevent certification of
an inferential contrast. Peak results remain accepted; no differential sex
effect is claimed. See the [result index](../README.md) for both arms.
