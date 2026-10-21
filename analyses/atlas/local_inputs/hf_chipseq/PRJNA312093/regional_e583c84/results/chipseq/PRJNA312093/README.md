# PRJNA312093 — ChIP-seq results

Peak processing is complete for the **H3K27me3 broad** and **H3K4me3 narrow**
arms. RNA-seq results and their admission decisions are independent and
unchanged. Differential Binding is **on HOLD**, not executed: see the
[comparison-specific eligibility review](../../../docs/chipseq/prjna312093_differential_binding.md).

The same preserved per-library MACS3 calls were consolidated using the
nf-core/chipseq `2.1.0` policy, BEDTools `2.30.0`, and **at least one calling
library**. HOMER `4.11` annotated the merged regional catalog against
SM_V10/WBPS19. There are two deposited IP libraries per group, with explicit
matched Input IDs recorded in each group manifest.

## H3K27me3

| Group | Regional candidates | Covered bases | Associated gene IDs |
|---|---:|---:|---:|
| Female adult | 75,155 | 19,954,038 | 6,120 |
| Female cercaria | 29,089 | 11,484,061 | 4,326 |
| Male adult | 39,845 | 9,901,810 | 5,031 |
| Male cercaria | 44,179 | 30,338,539 | 5,455 |

- [Machine-readable group summary](broad/regional_summary.tsv).
- Group catalogs and support: [female adult](broad/chipseq/consensus/female_adult/manifest.json),
  [female cercaria](broad/chipseq/consensus/female_cercaria/manifest.json),
  [male adult](broad/chipseq/consensus/male_adult/manifest.json), and
  [male cercaria](broad/chipseq/consensus/male_cercaria/manifest.json).
- Gene associations: [female adult](broad/chipseq/peak_annotation/female_adult/peak_gene_associations.tsv.gz),
  [female cercaria](broad/chipseq/peak_annotation/female_cercaria/peak_gene_associations.tsv.gz),
  [male adult](broad/chipseq/peak_annotation/male_adult/peak_gene_associations.tsv.gz), and
  [male cercaria](broad/chipseq/peak_annotation/male_cercaria/peak_gene_associations.tsv.gz).
- [Publication and comparison record](manifests/regional_publication.json)
  and [execution review](reports/PRJNA312093_execution_report.md).

Each group includes the unfiltered/selected regional BED, complete support
matrix, source-peak mapping, upstream boolean matrix, separate ≥2-library core,
normalized annotated peaks, gene associations and original HOMER rows. Tables
are gzip-compressed; hashes identify both the published files and native source
content. BED coordinates are zero-based and half-open.

## Evidence and limits

The minimum-one catalog retains regional candidates, including calls found
only in one library. Coverage of the input union is a consequence of this
policy, not measured biological recall. `single_library` and `all_libraries`
describe membership within these two-library groups. Coverage fractions and
the shared core remain separate; neither membership nor IP/Input ratios are
converted into statistical significance or verified reproducible occupancy.

`occupancy_status=not_assessed` is explicit. This source execution also predates
the native HOMER TSS snapshot, so gene-assignment ambiguity is `not_assessed`,
not `false`. Independent direct-tool comparison disagreed on one nearest-gene
assignment in each male group, involving `Smp_095580` and `Smp_341610`, which
share a TSS. The original strict `FAIL_DIRECT_EQUIVALENCE` record is retained;
the assignments were not changed to force agreement. Both female groups passed
the direct comparison. Nearest-gene associations do not establish regulation.

No Differential Binding is claimed. This compact catalog publication is not
a new complete native terminal run, does not certify all upstream QC, and is
not a self-contained re-entry package. FASTQs, BAMs, bigWigs, references,
workdirs and private execution logs are excluded. The source commit and exact
reference/annotation hashes are recorded per group; these runs must not be
relabeled as executions of a later HelixForge release.

## H3K4me3

| Group | Regional candidates (minimum one library) | Associated gene IDs | Regions retained by IDR ≤0.05 |
|---|---:|---:|---:|
| Female adult | 12,652 | 7,031 | 4,348 |
| Female cercaria | 12,946 | 6,933 | 4,461 |
| Male adult | 14,088 | 7,128 | 4,719 |
| Male cercaria | 14,151 | 7,103 | 4,651 |

- [Group summary](narrow/regional_summary.tsv),
  [per-library Peak QC](narrow/qc/peak_qc_summary.tsv.gz) and
  [coverage-track inventory](narrow/tracks/tracks.tsv.gz).
- Selected catalogs and gene associations: [female adult](narrow/chipseq/consensus/female_adult/manifest.json),
  [female cercaria](narrow/chipseq/consensus/female_cercaria/manifest.json),
  [male adult](narrow/chipseq/consensus/male_adult/manifest.json),
  [male cercaria](narrow/chipseq/consensus/male_cercaria/manifest.json).
- [Native continuation report](reports/PRJNA312093_H3K4me3.html),
  [MultiQC](reports/PRJNA312093_H3K4me3_multiqc.html) and
  [execution summary](narrow/manifests/native_execution_summary.json).

H3K4me3 was processed with HelixForge commit
`2833442eb01edbe2e090870ba945fc074bc4d564`, which includes the continuation
fix after v1.2.5. Four HOMER groups and 20 coverage tracks were produced;
the native terminal manifest passed schema, semantic and tracked-file
validation. The continuation report covers its supplied producers; upstream
alignment manifests were not included in that report inventory. This is not
a combined native report for both marks.

The **minimum-one regional catalog**, not the IDR-filtered set or shared core,
is the annotation input. IDR remains a separately identified replicate
reproducibility result at 0.05, ranked by signal value. It does not certify all
regions in the regional catalog. HOMER assignments, library memberships,
coverage fractions and shared core are preserved without reassignment.

Unlike the earlier H3K27me3 execution, the H3K4me3 source includes HOMER TSS
snapshots and native per-region ambiguity fields. `true`, `false` and
`not_assessed` states, with their original reasons, are retained; an unmatched
provider position/distance is not converted into an unambiguous assignment.
Neither a nearest-gene assignment nor a peak in one library proves regulation
or reproducible occupancy. Large BAMs/bigWigs and private audit data are not
distributed in Git; the public package is not a standalone re-entry bundle.

## Sources

[BioProject PRJNA312093](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA312093) and
[nf-core/chipseq 2.1.0](https://github.com/nf-core/chipseq/tree/2.1.0).
