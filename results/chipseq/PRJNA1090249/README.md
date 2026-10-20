# PRJNA1090249 — HP1 ChIP-seq results

Six HP1 IP libraries and two Inputs from cercariae and sporocysts of
*Schistosoma mansoni* were processed against SM_V10/WBPS19. Read processing
and per-library MACS3 `3.0.4` calls used HelixForge `v1.2.0`. The active
regional catalog uses the nf-core/chipseq `2.1.0` consolidation policy,
BEDTools `2.30.0`, support from **at least one library**, and HOMER `4.11`.

## Current results

- [Regional and gene-association summary](broad/regional_summary.tsv).
- [Cercarial catalog and evidence](broad/chipseq/consensus/cercaria/manifest.json)
  and [sporocyst catalog and evidence](broad/chipseq/consensus/sporocyst/manifest.json).
- [Cercarial peak–gene associations](broad/chipseq/peak_annotation/cercaria/peak_gene_associations.tsv.gz)
  and [sporocyst peak–gene associations](broad/chipseq/peak_annotation/sporocyst/peak_gene_associations.tsv.gz).
- [Per-library QC](qc/qc_sample_summary.tsv), [QC summary](qc/qc_summary.json),
  [MultiQC](reports/PRJNA1090249_multiqc.html), and unchanged per-library
  MACS3 calls under `broad/080-peak-calling/`.
- [Publication and comparison record](manifests/regional_publication.json),
  [execution review](reports/PRJNA1090249_execution_report.md), and
  [reference manifest](manifests/reference_bundle.manifest.json).

| Group | Regional candidates | Covered bases | Associated gene IDs |
|---|---:|---:|---:|
| Cercariae | 1,816 | 549,790 | 1,480 |
| Sporocysts | 867 | 268,249 | 771 |

These catalogs replace the former 367/98 atomic support segments as the
active input to gene association. Full merged regional intervals are retained;
the core with support from at least two libraries remains a separate spatial
product and does not crop the catalog. BED coordinates are zero-based and
half-open. Gene counts describe nearest-gene associations, not proven targets.

## Evidence and interpretation

The TSVs identify calling libraries and samples, their counts, support class,
per-library covered bases/fractions, and shared core. `single_library`,
`multiple_libraries`, and `all_libraries` describe library membership, not
biological reproducibility or statistical significance. All regions retain
`occupancy_status=not_assessed`.

These completed source executions predate the native HOMER TSS snapshot.
Their ambiguity fields explicitly say `not_assessed`; they must not be
interpreted as unambiguous. Gene assignments are unchanged. The independent
direct-tool comparison passed for both HP1 groups.

The 200 bp model-free extension remains an operational assumption, not a
measured fragment length. [Fragment diagnostics](qc/fragment_length_diagnostic/index.html)
and [article-inspired comparisons](qc/article_method_comparison/relatorio.html)
remain available. Stage/strain confounding and library-source independence
still limit inference; no Differential Binding was performed. The study
remains accepted for atlas inclusion, with evidence strength attached to each
association. See the [peak policy](../../../docs/chipseq/prjna1090249_peak_policy.md)
and [DB review](../../../docs/chipseq/differential_binding_review.md).

Large TSV/BED files are gzip-compressed. Original source hashes and public-file
hashes are recorded in the group manifests and `SHA256SUMS`. FASTQs, BAMs,
bigWigs, caches and private logs are excluded; this is not a self-contained
re-entry bundle. The superseded consensus, annotation and report were removed
from the active tree and remain recoverable through Git history.

## Sources

[BioProject PRJNA1090249](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1090249),
[study publication](https://doi.org/10.1590/0074-02760240075), and
[nf-core/chipseq 2.1.0](https://github.com/nf-core/chipseq/tree/2.1.0).
