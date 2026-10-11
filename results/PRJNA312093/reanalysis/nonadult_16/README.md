# PRJNA312093 non-adult DESeq2 refit

This supplementary analysis uses the validated gene-count matrix and sample
table from the original 20-library run. It includes 16 non-adult samples (two
female and two male biological replicates at each of four stages) and excludes
all four adult samples **only from this refit**. No FASTQ, Salmon quantification,
import artifact, original contrast, or terminal manifest was replaced.

The scientific settings remain `~ condition`, Wald test, rounded imported
counts, total-count filter `>=10`, `alpha=0.05`, and reported significance at
`padj < 0.05` with `|log2FC| >= 1`. DESeq2 version: `1.42.1`. The single-end
Salmon fragment prior for the source quantification is `200/80`; the separate
`250/80` sensitivity check passed but is not substituted here.

| Public artifact | Meaning |
| --- | --- |
| `contrasts/condition__*.tsv` | Four non-adult female-versus-male DESeq2 tables; these supersede the corresponding original tables for inference. |
| `included_samples.tsv`, `excluded_adults.tsv` | Exact sample scope; excluded adult libraries are retained in the original audit. |
| `nonadult_replicate_correlations.tsv`, `vst_correlation_matrix.tsv`, `nonadult_pca.tsv` | Exploratory VST concordance evidence. |
| `comparison_to_original.tsv` | Change in significant-gene sets and effect estimates relative to the original 20-sample model. |
| `decision.json` | Machine-readable acceptance and atlas selection contract. |

Eight non-adult replicate-pair VST correlations range from 0.9446 to 0.9937.
Across the four contrasts, 91.84–99.57% of genes significant in the original
fit remain significant, with no direction reversal among genes significant in
both fits. Significant-gene set membership is nevertheless sensitive to the
adult-sample exclusion, particularly at schistosomulum stages 1 and 3
(Jaccard 0.381 and 0.422). Thus the non-adult evidence is
`PASS_WITH_LIMITATIONS`; the adult comparison remains `NOT_EVALUABLE`.

The original 20-sample tables and terminal manifest remain available for
audit. A downstream atlas must select only the four tables here as inferential
Picard evidence and must not expose the original adult contrast as accepted.
The atlas has not yet been updated to enforce that selection.
