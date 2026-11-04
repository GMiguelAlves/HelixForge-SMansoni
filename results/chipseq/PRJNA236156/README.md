# PRJNA236156 — four-mark ChIP-seq results

Accepted processing-method reanalysis of **48 IP libraries and 12 matched
unbound-chromatin controls** from [Roquis et al. (2016)](https://doi.org/10.1111/mec.13555),
*Frequency and mitotic heritability of epimutations in Schistosoma mansoni*.
The panel covers H3K4me3, H3K27ac, H3K27me3 and H4K20me1 in cercariae and
adults, with sympatric/allopatric host contexts and the original Fa/Fb units.
The reference is **SM_V10/WBPS19**.

Standard HelixForge QC and Bowtie2 alignments precede the article-guided
SAMtools selection, Pyicos duplicate cap, subsampling and **PeakRanger 1.18
`ranger`** calls. The [analysis policy](../../../docs/chipseq/prjna236156_analysis.md)
documents the settings and adaptations. This is a processing reanalysis, not
an execution of the paper's Differential Binding or chromatin-state analyses.

## Final results

| Mark | IP libraries | Regions per library | Sum across libraries |
|---|---:|---:|---:|
| H3K4me3 | 12 | 7,752–14,808 | 125,089 |
| H3K27ac | 12 | 770–83,441 | 127,913 |
| H3K27me3 | 12 | 2,692–16,789 | 131,989 |
| H4K20me1 | 12 | 3,274–61,621 | 173,368 |

The **558,359** records are summed per-library calls, not unique genomic loci
or gene associations. No shared-core intersection or replicate consensus
recuts these intervals. All exported regions have the caller's `fdrPassed`
label. Individual library identity remains available for later interpretation.

- [Per-library peaks](peaks/peakranger/), preserving the called boundaries.
- [Sample/control table](metadata/sample_table.tsv) and
  [per-library QC](qc/qc_sample_summary.tsv).
- [Peak counts, coverage and IP/Input read fractions](metrics/peak_summary.tsv).
- [Comparison with Appendix S3](metrics/paper_comparison.tsv), matched by
  stage, host strain, clone, biological unit and mark.
- [QC summary](qc/qc_summary.json), [parameters](manifests/parameters.json),
  [versions](manifests/software_versions.tsv),
  [source checksums](manifests/source_peaksets.tsv), and
  [file inventory](manifests/files.json).

## Depth and interpretation

H3K4me3 cercaria calls use 16,000,000 selected reads per IP/Input where
available. Gua Fa1 uses all 14,730,063 eligible reads; Gua Fb5 uses all
5,925,466, close to the article's explicit 6,000,000 exception. H3K4me3 adult
calls retain the accepted 5,925,466-read depth. Other marks generally use
17,774,492 reads, with lower available depths for cercaria Gua Fb8 H3K27ac
and H4K20me1. Every actual depth is recorded per library.

For H3K4me3, mean region counts are 12,615 in cercariae versus 11,967 in the
article (+5.4%), and 8,233.2 in adults versus 8,377.5 (−1.7%). Similar counts
do not establish identical genomic positions. The other marks retain larger
differences from the published counts; no thresholds were changed to match
them.

Cercaria Gua Fb8 H3K27ac and H4K20me1 produce 83,441 and 61,621 regions,
respectively. These are the same libraries flagged as low-quality replicates
in Appendix S2. Their labels and QC metrics are retained, not silently pooled
or excluded. IP/Input fractions are descriptive QC, not new significance
tests or proof of reproducible occupancy.

Public files are compact BED3 gzip views of the native PeakRanger outputs,
with coordinates retained exactly and no liftover or offset transformation.
Native scores, summits, details, QC reports and processing provenance are
preserved separately. No new gene annotation, DB, IDR, EpiCSeg state assignment
or atlas update is included. FASTQs, BAMs, workdirs, execution scripts and
private paths/logs are excluded from Git; this is not a self-contained re-entry
package.
