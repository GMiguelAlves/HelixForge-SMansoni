# E-MTAB-451 provenance

Public, compact provenance for the accepted HelixForge RNA-seq execution.
Absolute server paths, raw FASTQs, the Nextflow work directory, and the private
audit archive are intentionally excluded from this repository.

- `execution_state.json`: frozen software, design, and execution summary.
- `final_validation.json`: scientific and operational acceptance gates.
- `cleanup.json`: post-acceptance retention and cleanup record.
- `performance.tsv`: sanitized Nextflow trace preserved with the result bundle.
- `metadata/E-MTAB-451/runs.tsv` and `samples.tsv`: the frozen 12-run to
  11-sample mapping.
- `source_snapshots/`: frozen public metadata used to curate the study.

The versioned results under `results/E-MTAB-451/` have their own
`SHA256SUMS`. The larger private work audit is retained outside Git and can be
used for controlled review if necessary.
