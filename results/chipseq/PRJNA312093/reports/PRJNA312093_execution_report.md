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

The H3K4me3 arm, whole-study terminal report, final combined QC review and DB
eligibility review are not certified by this publication. See the
[result index](../README.md) and group manifests for provenance and artifact hashes.
