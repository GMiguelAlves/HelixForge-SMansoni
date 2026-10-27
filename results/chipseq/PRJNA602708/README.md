# PRJNA602708 — ChIP-seq results

Separate `mapq30` and `mapq0` reanalysis layers; neither is certified occupancy or an exact reproduction of the article. Individual Inputs, per-library duplicate removal, mitochondrial read/mate exclusion and MACS3 p=0.01 are held fixed. The weak broad threshold is p=0.1, not FDR.

Each layer/mark contains native per-library peaks, nf-core/chipseq2.1.0 regional catalogs (minimum one library), separate shared core, HOMER4.11 annotations with library support/ambiguity, peak QC and a native HTML report. These mapping layers must not be merged without identifying their policies.

BAMs/BAIs, BigWigs and complete private audit are preserved separately, not committed. This is not a self-contained re-entry package. Differential Binding is suspended; no new IDR result is claimed.

Read [the analysis policy](../../../docs/chipseq/prjna602708_analysis.md) before interpreting the [regional summary](regional_summary.tsv) or [QC summary](qc/qc_summary.json).
