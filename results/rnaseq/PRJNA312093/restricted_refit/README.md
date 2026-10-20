# PRJNA312093 — approved restricted refit for atlas use

This is the reviewed **12-library** DESeq2 refit, not the frozen 20-library
workflow's primary DE output. It contains two biological replicates for each
female/male group at schistosomulum stages 1, 2 and 3. The original
`~ condition` model, total-count filter (`>=10`), Wald test and reporting
criteria (`padj < 0.05`, `|log2FC| >= 1`) were retained. The six within-group
VST replicate correlations range from 0.9415 upward.

Only these three tables in [`contrasts/`](contrasts/) are approved as atlas
DE sources. They contain 9,342 tested genes each and 290, 1,045 and 437
significant genes respectively. Their SHA-256 values match the preserved
private refit audit exactly; [`manifest.json`](manifest.json) records the
specific source files and eligibility warnings without private paths.

The original five contrast tables in `../differential_expression/contrasts/`
remain technical execution evidence, **not** atlas inputs. Adult and cercaria
sex contrasts remain `NOT_EVALUABLE` because of unresolved source-library
profile inconsistencies. No library was relabelled. The approved tables
must be displayed with `RESTRICTED_REFIT`, `TWO_REPLICATES_PER_GROUP` and
`SOURCE_IDENTITY_LIMITATION` warnings. This admission does not make the
20-sample analysis biologically accepted. If showing this study's TPM in the
atlas, select only the 12 columns named in `included_samples.tsv` from the
original expression matrix; do not silently display the eight adult/cercaria
libraries beside the restricted contrasts.

The [atlas admission record](../../../../docs/rnaseq/atlas_study_decisions.md)
contains the investigation and comparison to the intermediate 16-sample
sensitivity fit. The sample list and replicate correlations are included
here; R session details and scheduler logs remain in private audit storage.
