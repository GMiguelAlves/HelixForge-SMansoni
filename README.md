# HelixForge-SMansoni

Reproducible RNA-seq, ChIP-seq, and cross-omic integration analyses of
*Schistosoma mansoni* using
[HelixForge](https://github.com/GMiguelAlves/HelixForge).

## Explore the atlas

The [atlas release](https://github.com/GMiguelAlves/HelixForge-SMansoni/releases/tag/atlas-v2026.10.21) provides a ready-to-open RNA/ChIP atlas. Extract `atlas-init.zip` and run `python atlas.py init` for Docker, or extract `atlas-local.zip` and run `python atlas.py init --local` without Docker. Both open **http://127.0.0.1:8765/**. Python 3.10+ is required; no analysis dependencies are needed.

From a clone, the same commands work directly; `docker compose up --build -d` also builds and starts the atlas. See [installation, scientific scope and size notes](ATLAS.md). The frozen atlas is a consultation layer over the accepted results below; it does not rerun or broaden their scientific analyses.

```text
PROJECT_STATUS = ACTIVE
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

Study- and assay-specific evidence enters shared analyses only after its
scientific design, metadata and execution contracts are reviewed.

| Runtime baseline | Version |
|---|---|
| Nextflow | `25.10.7` |
| Java | `21` |

Exact HelixForge releases and commits are recorded per study, not imposed as
one project-wide pin; see the study registry and execution provenance.

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
The pinned HelixForge runtime validates and consumes that immutable index
explicitly; it does not depend on `-resume` to avoid rebuilding it.

## Studies and preparation

Study status, execution reviews and ChIP-seq preparation are documented in the
[study registry](docs/project/studies_and_preparation.md).

## Execution policy

Only one study may occupy heavy scratch storage at a time. Data acquisition is
outside the scientific Nextflow DAG, and all heavy work on an HPC system must
run inside Slurm allocations. Persistent results, terminal manifests,
provenance, and compact audit packages are retained; FASTQs, work directories,
indexes, and other reproducible intermediates are removed after acceptance.

See [server execution](docs/operations/server_execution.md) and the
[server preparation report](docs/operations/server_preparation_report.md) before staging
any data.

## Repository structure

- `config/rnaseq/` and `config/chipseq/`: assay-specific study settings;
  shared reference and site configuration stay at `config/`.
- `metadata/rnaseq/` and `metadata/chipseq/`: assay-specific curated metadata,
  study inventories and input contracts.
- `scripts/`: download, validation, conservative harmonization, provenance,
  and study launch helpers.
- `results/rnaseq/` and `results/chipseq/`: compact, versionable assay-specific
  results, QC, reports and sanitized provenance. No raw or heavy intermediates.
- `provenance/rnaseq/`: per-study RNA-seq records; shared reference and frozen
  source-package evidence stay at `provenance/`.
- `analyses/`: cross-study, RNA-seq/ChIP-seq integration, atlas, and candidate
  analyses introduced under separately reviewed contracts.
- `docs/`: assay-specific, candidate, operational, and project documentation
  organized under [`rnaseq/`](docs/rnaseq/), [`chipseq/`](docs/chipseq/),
  [`candidates/`](docs/candidates/), [`operations/`](docs/operations/) and
  [`project/`](docs/project/).

## Reproducibility

The curated inputs were validated against the HelixForge contracts and retain
their original source checksums under `provenance/source_package/frozen/`. Each real
execution must record the HelixForge tag and commit, reference checksums,
commands, environment, Slurm metadata, terminal manifest, and an execution
summary under `provenance/rnaseq/<study>/` (or the corresponding assay folder).

Completion of an individual study does not automatically authorize a new
cross-study synthesis, coexpression analysis, candidate-prioritization pass, or
epigenomic integration. Each requires its own reviewed execution contract.

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
