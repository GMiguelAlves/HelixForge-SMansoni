# PRJEB3190 provenance

Public, compact provenance for the accepted HelixForge RNA-seq execution.
Absolute server paths, raw FASTQs, the Nextflow work directory, and the private
audit archive are intentionally excluded from this repository. The complete
candidate-gene report is also retained outside Git under the current
repository-size policy.

- `execution_state.json`: frozen software, design, and execution summary.
- `final_validation.json`: scientific and operational acceptance gates.
- `cleanup.json`: post-acceptance retention and cleanup record.
- `performance.tsv`: sanitized Nextflow trace preserved with the result bundle.
- `metadata/PRJEB3190/runs.tsv` and `samples.tsv`: the frozen 20-run to
  20-sample mapping.
- `source_snapshots/`: frozen public metadata used to curate the study.

The versioned results under `results/PRJEB3190/` have their own `SHA256SUMS`.
The larger private work audit is retained outside Git and can be used for
controlled review if necessary.
