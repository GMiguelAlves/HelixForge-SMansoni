# PRJNA1090249 — HP1 ChIP-seq results

Native ChIP-seq processing is complete for six HP1 IP libraries and two Inputs
from cercariae and sporocysts of *Schistosoma mansoni*. The analysis uses
HelixForge `v1.2.0`, SM_V10/WBPS19 and MACS3 `3.0.4` broad peak calling.
The study is accepted for inclusion in the atlas. Differential Binding is
not performed in the current analysis.

## Results

- [Standardized per-library QC](qc/qc_sample_summary.tsv) and
  [QC summary](qc/qc_summary.json).
- [Native ChIP-seq HTML report](broad/chipseq/report/report_result/chipseq_report.html)
  and [MultiQC report](reports/PRJNA1090249_multiqc.html).
- [Per-library peak QC](broad/chipseq/peak_qc/peak_qc_summary.tsv),
  [consensus summary](broad/chipseq/consensus/consolidation_summary.tsv),
  [annotated peaks](broad/chipseq/peak_annotation/peak_annotation_aggregate/annotated_peaks.tsv)
  and [peak–gene associations](broad/chipseq/peak_annotation/peak_annotation_aggregate/peak_gene_associations.tsv).
- [Fragment-length diagnostics and figures](qc/fragment_length_diagnostic/index.html)
  and [article-inspired comparison](qc/article_method_comparison/relatorio.html).
- [Final execution review](reports/PRJNA1090249_execution_report.md),
  [publication validation](manifests/publication_validation.json),
  [reference manifest](manifests/reference_bundle.manifest.json) and
  [execution provenance](../../../provenance/chipseq/PRJNA1090249/execution_provenance.json).

Consensus retained **367 cercarial and 98 sporocyst support segments** with
support from at least two deposited IP libraries. These are support segments,
not independently verified broad biological domains. Per-library peak counts
range from 48 to 1,240 and FRiP from 0.201% to 3.161%.

## Interpretation

The 200 bp model-free extension is an explicit operational assumption, not a
measured fragment length. The 150/200/250 bp sensitivity analysis and pooled
article-inspired follow-up are recorded in the
[peak-calling policy](../../../docs/chipseq/prjna1090249_peak_policy.md).
The follow-up did not establish exact reproduction of the publication.
TES-associated coverage was also present in Input; profile shape alone must
not be interpreted as HP1-specific enrichment.

Cercariae SmVEN and sporocysts SmDFO differ in both stage and strain. Source
independence of the three deposited libraries per group is not established.
These limitations remain attached to the atlas evidence; atlas inclusion is
not evidence of differential occupancy or causality. The
[final DB review](../../../docs/chipseq/differential_binding_review.md)
records why no inferential comparison was run.

## Public and retained artifacts

This reduced public package contains peaks, annotations, QC, figures, reports
and sanitized manifests. FASTQs, BAM/BAI, bigWigs, workdirs and private audit
logs are excluded from Git. Full persistent results and private audit are
retained separately. The native terminal manifest is retained as a provenance
record, but this reduced public copy is **not a self-contained re-entry
package**, because large referenced artifacts are omitted.

`TECHNICAL_EXECUTION = PASS`

`PROCESSING_STATUS = COMPLETE`

`ATLAS_INCLUSION = ACCEPTED`

`DIFFERENTIAL_BINDING = NOT_ELIGIBLE_IN_CURRENT_ANALYSIS`

## Source

[BioProject PRJNA1090249](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1090249)
and [study publication](https://doi.org/10.1590/0074-02760240075).
