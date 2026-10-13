# PRJEB3190

- Included paired-end runs: **20**.
- Biological samples: **20**.
- Samples requiring run/lane consolidation: **0**.
- Excluded ENA runs: **36**.

Only the canonical poly(A)/mRNA series is included: time 0 h has two
replicates; the other six time points have three each (20 samples total).
miRNA/small-RNA libraries and three noncanonical labels are excluded from
the `full_length` path and remain in the exclusion inventory.

`sequencing_batch_10363` and `sequencing_batch_10376` were derived from ENA
`experiment_alias` prefixes. ENA defines aliases as submitter-provided names;
the paper does not establish those prefixes as technical batches. They remain
available for PCA, diagnostics, and future investigation, but not as an
inferential covariate. The frozen design remains `~ condition` with this
warning.

Repeated `sample_id` values in `metadata.csv` are intentional. HelixForge
processes each run, consolidates FASTQs at sample level before quantification,
and uses the frozen contrasts in `de_spec.json` rather than inferring them
from names.
