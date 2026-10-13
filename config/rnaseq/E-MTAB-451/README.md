# E-MTAB-451

- Included paired-end runs: **12**.
- Biological samples: **11**.
- Samples requiring run/lane consolidation: **1**.
- Excluded ENA runs: **0**.

Inferential differential expression is restricted to the three replicated
groups. `adult_mixed` and `cercarial_tail` (n=1) remain available for
quantification and exploration but are absent from inferential contrasts.
`somule1` has two technical runs that must be consolidated.

Repeated `sample_id` values in `metadata.csv` are intentional. HelixForge
processes each run, consolidates FASTQs at sample level before quantification,
and uses the frozen contrasts in `de_spec.json` rather than inferring them
from names.
