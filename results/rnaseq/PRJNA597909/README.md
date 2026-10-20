# PRJNA597909 results

The first complete native RNA-seq application in this project is classified
`PASS_WITH_LIMITATIONS`. QC, Salmon 1.10.3 quantification, tximport under the
frozen `production_v1` policy, four DESeq2 contrasts, figures, persistent
publication, terminal-manifest validation, and cleanup passed. The limitations
are operational: terminal-only manifest recovery was required after the host
selected a Python without `jsonschema`, and a controlled `-resume` attempt did
not reuse eligible cache entries on shared NFS. No scientific stage was
recomputed during recovery.

- `qc/`: compact per-sample read retention, FastQC, and Salmon mapping summary.
- `expression/`: gene-level counts, TPM, effective lengths, import statistics,
  versions, and sanitized R session information.
- `differential_expression/`: normalized counts, four contrast tables,
  aggregate summaries, coefficients, dispersions, and scientific figures.
- `rnaseq/`: portable terminal manifest and its Integration API artifact
  bundle.
- `manifests/`: terminal validation record and a public-location manifest copy.
- `reports/`: final Markdown report, MultiQC, Nextflow execution/timeline/DAG
  HTML reports, scientific figures, performance and scratch summaries,
  cleanup record, and review gate. Per-read FastQC HTML files remain outside
  Git because their information is consolidated by MultiQC.

The authoritative portable manifest is
[`rnaseq/rnaseq_run_manifest.json`](rnaseq/rnaseq_run_manifest.json); its
relative artifact locations resolve entirely inside the versioned bundle.

`PRJNA597909_RNASEQ_ANALYSIS = PASS_WITH_LIMITATIONS`

`READY_FOR_PRJNA597909_REVIEW`
