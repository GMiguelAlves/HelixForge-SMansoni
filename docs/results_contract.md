# Public RNA-seq results contract

This contract separates three objects that previously became mixed across
historical studies:

1. the native HelixForge result tree, retained in persistent project storage;
2. the private audit/re-entry package, which may contain complete Salmon and
   Integration API artifacts;
3. the compact public result package versioned in this repository.

Only the third object is described here.

## Stable layout

```text
results/<study>/
├── README.md
├── qc/
│   ├── qc_sample_summary.tsv
│   └── qc_summary.json
├── expression/
│   ├── counts_matrix.tsv
│   ├── tpm_matrix.tsv
│   └── length_matrix.tsv
├── differential_expression/        # optional by scientific design
├── manifests/
│   ├── rnaseq_run_manifest.json
│   └── terminal_manifest_validation.json
├── reports/
└── gene_report/                    # optional presentation product
```

Additional files may be published inside these directories, but publication
must not create alternative locations for the stable artifacts. In particular,
MultiQC belongs in `reports/` when it is included in Git; the two compact QC
contract files always remain in `qc/`.

## Required interface

Every registered study must publish:

- the study README;
- both compact QC files;
- the three gene-level expression matrices;
- the terminal RNA-seq manifest and validation record;
- a results classification and any limitations in the README.

Differential-expression outputs are required only when the frozen design
supports inference. Execution reports, figures and MultiQC are strongly
recommended, but may remain private when their size or content is unsuitable
for Git; that decision must be recorded in the study README.

## Excluded from Git

The following are preserved outside the public package:

- FASTQs and per-read FastQC HTML/ZIP files;
- complete Salmon quantification directories;
- duplicated `integration_artifacts/` trees;
- workdirs, caches, references, indexes and scheduler logs;
- RDS objects and other reproducible heavy intermediates;
- reports deliberately excluded by the project publication policy.

A terminal manifest in the compact public package is provenance. Re-entry is
certified against the separately preserved audit package; Git is not used as a
substitute for that archive.

## Historical packages

Older accepted studies contain transitional layouts. They remain traceable in
Git history, but new publications must use this contract. Their compact QC
interface has already been normalized; the remaining physical-layout migration
should be performed as one dedicated maintenance change rather than silently in
an unrelated scientific PR.
