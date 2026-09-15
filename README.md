# HelixForge-SMansoni

Reproducible RNA-seq, ChIP-seq, and cross-omic integration analyses of
*Schistosoma mansoni* using
[HelixForge](https://github.com/GMiguelAlves/HelixForge).

```text
PROJECT_STATUS = PROJECT_INITIALIZED
DATA_ANALYSIS_STATUS = DATA_ANALYSIS_NOT_YET_COMPLETE
```

## Relationship to HelixForge

HelixForge is the organism-agnostic workflow software. This repository is a
scientific project that pins and uses a released HelixForge version; it is not
a specialized fork and does not copy the HelixForge core.

## Scientific scope

This is the species-specific analysis repository for the complete project. It
will hold the registered RNA-seq analyses, ChIP-seq analyses, and integration
of evidence across both modalities. Each dataset remains independently
configured and auditable, while accepted terminal manifests provide the
boundary for later integration.

The studies currently registered below are the first RNA-seq execution series.
Their presence does not restrict the repository to transcriptomics. ChIP-seq
datasets and integrative analyses are added only after their own scientific
design, metadata, and execution contracts are reviewed.

| Dependency | Frozen value |
|---|---|
| HelixForge release | [`v1.0.0`](https://github.com/GMiguelAlves/HelixForge/releases/tag/v1.0.0) |
| HelixForge commit | `14dc5a75d6c63f20d10c135f0c75138ea76dcc12` |
| Nextflow | `25.10.7` |
| Java | `21` |

## Reference

| Field | Frozen value |
|---|---|
| Organism | *Schistosoma mansoni* |
| Assembly | `SM_V10` |
| Source | WormBase ParaSite |
| Release | `WBPS19` (March 2024) |
| Reference ID | `Schistosoma_mansoni_SM_V10_WBPS19` |

Raw FASTA, transcriptome, GFF3, GTF, and Salmon indexes are not committed.
Their sources, filenames, release identity, retrieval date, and verified
checksums are frozen in
[`config/reference/reference_manifest.json`](config/reference/reference_manifest.json).
The validation metrics and reusable-index provenance are recorded under
[`provenance/reference/`](provenance/reference/).

## Initial studies

| Order | Study | Runs | Biological samples | Initial execution |
|---:|---|---:|---:|---|
| 1 | `PRJNA602528` | 10 | 10 | operational calibration through Import; DE blocked by design |
| 2 | `PRJNA597909` | 20 | 20 | full RNA-seq |
| 3 | `PRJEB14695` | 138 | 23 | full RNA-seq |
| 4 | `PRJEB32839` | 150 | 75 | full RNA-seq |

`PRJNA602528` is the first real *S. mansoni* application and operational
calibration run. It is not a new benchmark of HelixForge. Its public design has
one library per time point, so the frozen project stops after tximport and does
not perform DESeq2 inference.

## Execution policy

Only one study may occupy heavy scratch storage at a time. Data acquisition is
outside the scientific Nextflow DAG, and all heavy work on an HPC system must
run inside Slurm allocations. Persistent results, terminal manifests,
provenance, and compact audit packages are retained; FASTQs, work directories,
indexes, and other reproducible intermediates are removed after acceptance.

See [server execution](docs/server_execution.md) and the
[server preparation report](docs/server_preparation_report.md) before staging
any data.

## Repository structure

- `config/`: reference lock, site template, and study-specific analysis files.
- `metadata/`: run-level and sample-level curated metadata plus controlled
  vocabularies.
- `scripts/`: download, validation, conservative harmonization, provenance,
  and study launch helpers.
- `results/`: versionable summaries, tables, figures, reports, and terminal
  manifests; no raw or heavy intermediate data.
- `provenance/`: frozen input-package evidence and per-study execution records.
- `analyses/`: cross-study, RNA-seq/ChIP-seq integration, atlas, and candidate
  analyses introduced under separately reviewed contracts.
- `docs/`: methods, operational policy, dictionary, and limitations.

## Reproducibility

The curated inputs were validated against the HelixForge contracts and retain
their original source checksums under `provenance/source_package/frozen/`. Each real
execution must record the HelixForge tag and commit, reference checksums,
commands, environment, Slurm metadata, terminal manifest, and an execution
summary under `provenance/<study>/`.

This initialization does not start downloads, processing, meta-analysis,
atlas construction, coexpression, candidate prioritization, or epigenomic
integration.

## Administrative validation

The repository checks metadata cardinality and identity, the conservative
master manifest, reference identity, the frozen source-package checksums,
shell syntax, JSON syntax, and exclusion of sequencing-scale files. Run:

```bash
python -m unittest discover -s tests -p 'test_*.py'
python scripts/validate/validate_project.py
```

## Citation and license

Project citation metadata is in [`CITATION.cff`](CITATION.cff). Code,
configuration, and documentation are licensed under Apache License 2.0. Public
datasets and reference resources retain their original terms and attribution.
