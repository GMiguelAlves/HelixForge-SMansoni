# PRJNA312093 — H3K27me3 regional results

This publication contains the completed **H3K27me3 broad arm**, not the entire
ChIP-seq study. The H3K4me3 arm is pending. RNA-seq results and their admission
decisions are independent and unchanged.

The same preserved per-library MACS3 calls were consolidated using the
nf-core/chipseq `2.1.0` policy, BEDTools `2.30.0`, and **at least one calling
library**. HOMER `4.11` annotated the merged regional catalog against
SM_V10/WBPS19. There are two deposited IP libraries per group, with explicit
matched Input IDs recorded in each group manifest.

## Results

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

## Sources

[BioProject PRJNA312093](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA312093) and
[nf-core/chipseq 2.1.0](https://github.com/nf-core/chipseq/tree/2.1.0).
