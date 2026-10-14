# PRJNA428470 — RNA-seq execution and interpretation review

## Design and execution

Six paired-end runs represent six independent adult-female biological
libraries, three exposed to 5-azacitidine and three controls. The frozen
[DE specification](../../../../config/rnaseq/PRJNA428470/de_spec.json) uses
`~ condition`, a Wald test, total count at least 10, `padj < 0.05` and
`|log2FC| >= 1`. No sample, contrast, scientific parameter or threshold was
changed during review. HelixForge `v1.1.0`, Nextflow `25.10.7`, Salmon
`1.10.3` and reference SM_V10/WBPS19 completed 77/77 tasks. All 12 accepted
FASTQs passed ENA size/MD5, SHA-256, gzip and paired-read checks; two
corrupt copies from an initial transfer were discarded **before** analysis
and replaced by verified downloads.

## Quality-control gate

| Measure | Result |
| --- | ---: |
| Biological libraries | 6 (3 treated, 3 control) |
| Post-trim read retention | 98.62%–99.18% |
| Salmon mapping to mRNA index | 34.61%–60.83% |
| Libraries below 60% mapping review threshold | 4/6 |
| Mean mapping, treated / control | 47.08% / 50.11% |
| Minimum within-group log2-CPM correlation | 0.9888 |
| MultiQC / terminal manifest | `PASS` / validated |

The [QC summary](../qc/qc_summary.json), [per-library table](../qc/qc_sample_summary.tsv)
and [replicate diagnostics](../qc/replicate_qc.tsv) preserve these
classifications. FastQC base and per-sequence quality passed for all raw
mates, while duplication and base-composition modules raised warnings/fails;
those flags alone were not used to remove RNA-seq libraries. Mapping is a
real scope limitation, not an indication that the computational workflow
failed.

## Mapping investigation and biological check

A deterministic approximately 0.5% paired-read subset from one
higher-mapping and one lower-mapping treated run reproduced their Salmon
mapping difference. Independent STAR alignment to the SM_V10 genome found
fewer unique and more highly multimapping reads in the lower-mapping run;
a limited ribosomal-reference screen also found more rRNA-compatible R1
reads. This supports a library-composition contribution in **two** runs,
not an exact six-library rRNA percentage or proof of cross-species
contamination. The original [Amaral et al. (2020) reanalysis](https://www.nature.com/articles/s41598-020-78669-5)
included both coding transcripts and lncRNAs under v7/WBPS14. Our Salmon
index contains SM_V10/WBPS19 mRNAs, so absolute mapping rates and lncRNA
results are not directly comparable.

The article reported upregulation of `Smp_151640` by approximately 14.6-fold
and downregulation of `Smp_121390` by approximately 3.6-fold in the same
RNA-seq dataset. Our coding-gene contrast has the same directions and
approximately 14.9-fold up / 3.6-fold down. This is a targeted plausibility
check, **not** independent validation of every significant gene.

## Differential-expression stability

The full, frozen six-library DESeq2 fit tested 9,307 genes and marked 1,454
significant. A supplementary sensitivity analysis omitted each library in
turn, refitted the same model with the same count filter, and compared each
five-library result to the full fit:

| Diagnostic | Range across six omissions |
| --- | ---: |
| Spearman correlation of gene log2FC | 0.977–0.991 |
| Originally significant genes remaining significant | 1,286–1,363 / 1,454 (88.4%–93.7%) |
| Direction agreement among originally significant genes present in each refit | 100% |

The [per-omission table](../qc/leave_one_out_summary.tsv) is an exploratory
robustness check; it is not a new primary model, and no library was excluded
from the published [contrast table](../differential_expression/contrasts/DEG_condition__female_adult_5azac_vs_female_adult_control.tsv).

## Decision and provenance

`TECHNICAL_EXECUTION=PASS`, `SAMPLE_QC=REVIEW_4_OF_6` and
`DE_CODING_GENE_SCOPE=PASS_WITH_LIMITATIONS`. The frozen contrast may be
displayed in the atlas only with a prominent mapping warning and only for
coding genes represented by the declared mRNA index. The 1,454 genes are
statistical candidates, not confirmed treatment targets. No inference is
made about lncRNAs, unquantified reads or between-study batch effects.

The original result tree and independent mapping investigation were
archived privately before study-specific scratch cleanup. A second verified
private ZIP contains the six leave-one-out fits, R session information,
script, scheduler record and dated Portuguese decision. The public
[terminal manifest](../manifests/rnaseq_run_manifest.json) and
[validation record](../manifests/terminal_manifest_validation.json) refer
to the complete original tree; unavailable artifacts in this compact
package remain identified by provenance rather than silently replaced.
