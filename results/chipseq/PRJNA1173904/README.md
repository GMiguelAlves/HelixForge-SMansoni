# PRJNA1173904 — eight-mark gonad ChIP-seq results

The [full panel](full_panel/README.md) contains **64 IP libraries**, eight marks,
four gonad/mating-experience conditions and two labelled replicates per group,
processed against **SM_V10/WBPS19**. All libraries are included **with warnings**;
this acceptance does not certify every called region as true occupancy.

- [Per-library ChromstaR calls](full_panel/peaks/) retain the complete native
  interval boundaries; no replicate consensus recuts the regions.
- [HOMER gene annotations](full_panel/annotation/) retain library, condition,
  mark, control policy, QC/review warnings and gene-assignment ambiguity.
- [Samples](full_panel/metadata/sample_table.tsv),
  [QC](full_panel/qc/qc_sample_summary.tsv),
  [model metrics](full_panel/metrics/model_summary.tsv),
  [annotation metrics](full_panel/metrics/annotation_summary.tsv),
  [replicate overlap](full_panel/metrics/replicate_region_overlap.tsv),
  [fit figures](full_panel/figures/) and
  [parameters](full_panel/manifests/parameters.json) describe the full panel.

There are **1,923,296 per-library regional records**, not unique genomic loci
or validated regulatory interactions. Ovary IPs use the Input for their mating
condition. Testis IPs are a **separate no-Input alternative** because the two
testis Inputs are restricted; they are not a controlled reproduction of the
article. Discordant H3K4me3/H3K9me3 testis libraries remain explicitly flagged.
See the [processing policy](../../../docs/chipseq/prjna1173904_analysis.md).

## Earlier H3K4me3 ovary calibration panel

The results below describe the earlier four-IP calibration, retained separately.
Its downsampling, common Input and three-tool comparison do not describe the
new full panel.

Accepted processing-method panel following the H3K4me3 ovary procedure in
[Augusto et al. (2025)](https://doi.org/10.3390/ncrna11020025), adapted to
**SM_V10/WBPS19**. It covers **four H3K4me3 IP libraries**, two per mating
experience, with one common Input. It is not an eight-mark reanalysis or a
chromatin-state atlas of the complete BioProject.

The study's published **methods**, rather than its downstream biological
analyses or target peak counts, guide this processing. Existing QC and
Bowtie2 alignments were reused. Each deduplicated IP contributes 2,500,000
complete fragments; the common Input contributes all 4,509,411 eligible
fragments. Reference, sampling-unit, software and control-policy adaptations
are explicit in the [analysis policy](../../../docs/chipseq/prjna1173904_analysis.md).

### Calibration results

| Library | Mating experience | Replicate | MACS2 regions | Peakranger regions | ChromstaR univariate regions |
|---|---|---:|---:|---:|---:|
| SRR31021582 | Experienced | 2 | 8,064 | 6,901 | 6,779 |
| SRR31021583 | Experienced | 1 | 7,492 | 9,607 | 7,555 |
| SRR31021601 | Inexperienced | 2 | 8,136 | 6,459 | 7,409 |
| SRR31021602 | Inexperienced | 1 | 7,743 | 5,579 | 8,034 |

- Per-library intervals: [MACS2](peaks/macs2/),
  [Peakranger](peaks/peakranger/), and
  [ChromstaR univariate](peaks/chromstar_univariate/).
- [Sample/control table](metadata/sample_table.tsv),
  [per-library QC](qc/qc_sample_summary.tsv), and
  [review summary](qc/qc_summary.json).
- [MACS2/Peakranger counts and coverage](metrics/caller_peak_counts.tsv),
  [their interval concordance](metrics/caller_overlap.tsv), and
  [ChromstaR counts/convergence](metrics/chromstar_univariate_summary.tsv).
- [ChromstaR interval comparisons](metrics/chromstar_overlap.tsv) and
  [explicit chromosome-name mapping](metrics/chromstar_chromosome_mapping.tsv).
- [Processing parameters](manifests/parameters.json),
  [software versions](manifests/software_versions.tsv),
  [source peakset checksums](manifests/source_peaksets.tsv), and
  [file inventory](manifests/files.json).

All distributed peak files are **BED3, zero-based, half-open**, compressed with
gzip. They preserve individual called interval boundaries; they are not merged
consensus sets. ChromstaR export aliases were mapped explicitly to SM_V10
sequence names without shifting coordinates. Native scores, complete native
outputs and fitted objects are preserved separately, not reconstructed from
these BED3 files.

### Calibration replicate concordance

Jaccard measures overlap of covered bases, not biological accuracy.

| Ovary group | MACS2 | Peakranger | ChromstaR univariate |
|---|---:|---:|---:|
| Mating-experienced | 0.685864 | 0.579643 | 0.712139 |
| Mating-inexperienced | 0.549628 | 0.555583 | 0.621837 |

ChromstaR has higher replicate concordance in this panel, but also broader
coverage. Its regions contain 99.83–99.99% of MACS2-covered bases; additional
coverage is not independently established as true occupancy. Only separately
fitted **univariate** models were used for replicate comparisons, not the
jointly inferred per-replicate exports.

### Calibration fit diagnostics

All four univariate models and the joint model met their native convergence
criteria. Convergence alone does not establish biological validity.

| Library | Native distribution fit |
|---|---|
| SRR31021582 | [View figure](figures/SRR31021582_chromstar_fit.png) |
| SRR31021583 | [View figure](figures/SRR31021583_chromstar_fit.png) |
| SRR31021601 | [View figure](figures/SRR31021601_chromstar_fit.png) |
| SRR31021602 | [View figure](figures/SRR31021602_chromstar_fit.png) |

The panel is accepted as the method-following result. No new gene annotation,
Differential Binding, IDR or eight-mark combinatorial result is claimed.
Native final interval tables, fit figures, fitted ChromstaR objects and
processing metadata are retained separately. Earlier processing attempts and
their intermediate alignments are not part of the retained result set.
FASTQs, BAMs, bigWigs, execution scripts, private logs and workdirs are excluded
from Git. This is not a standalone re-entry package.
