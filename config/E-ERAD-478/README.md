# E-ERAD-478

- Included paired-end runs: **120**.
- Biological samples: **60**.
- Samples requiring run/lane consolidation: **60**.
- Excluded ENA runs: **20**.

The paper's corrected S11/S12 spreadsheet determines inclusion: 120 runs,
60 samples, and two runs per sample. The design adjusts for `batch` and tests
sex, pairing status, and time through explicit contrasts.

Repeated `sample_id` values in `metadata.csv` are intentional. HelixForge
processes each run, consolidates FASTQs at sample level before quantification,
and uses the frozen contrasts in `de_spec.json` rather than inferring them
from names.
