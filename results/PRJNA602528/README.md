# PRJNA602528 results

The first operational calibration run completed under the frozen
`IMPORT_ONLY` contract and is classified `PASS_WITH_LIMITATIONS`. The
limitations are operational: top-level `-resume` did not reuse eligible tasks
on shared NFS, and scratch was measured at checkpoints rather than
continuously. QC, Salmon quantification, tximport, persistent publication, and
terminal-manifest validation passed.

- `qc/`: compact per-sample retention and Salmon mapping summary.
- `expression/`: gene-level counts, TPM, effective-length matrix, and import
  summaries.
- `compatibility/quants/`: the ten compact `quant.sf` files required by the
  portable terminal manifest.
- `manifests/`: validated Integration API terminal manifest.
- `reports/`: execution report and descriptive, non-inferential figures.

`DIFFERENTIAL_EXPRESSION = NOT_PERFORMED_BY_DESIGN`

The public replication structure does not support the frozen differential-
expression design. No contrast, DEG table, p-value, or FDR claim was produced.
