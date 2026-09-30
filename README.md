# HelixForge-SMansoni

Reproducible RNA-seq, ChIP-seq, and cross-omic integration analyses of
*Schistosoma mansoni* using
[HelixForge](https://github.com/GMiguelAlves/HelixForge).

```text
PROJECT_STATUS = ACTIVE
DATA_ANALYSIS_STATUS = RNASEQ_ATLAS_READY_FOR_REVIEW
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
| HelixForge release | [`v1.0.2`](https://github.com/GMiguelAlves/HelixForge/releases/tag/v1.0.2) |
| HelixForge commit | `5d4b3e696319db5cd7633472504964f1dc7c0434` |
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
HelixForge `v1.0.2` validates and consumes that immutable index explicitly;
it does not depend on `-resume` to avoid rebuilding it.

## Initial studies

| Order | Study | Runs | Biological samples | Initial execution |
|---:|---|---:|---:|---|
| 1 | `PRJNA602528` | 10 | 10 | Import completed; `PASS_WITH_LIMITATIONS`; DE not applicable by design |
| 2 | `PRJNA597909` | 20 | 20 | Full RNA-seq completed; `PASS_WITH_LIMITATIONS`; DE enabled |
| 3 | `PRJEB14695` | 138 | 23 | Full RNA-seq completed; `PASS_WITH_LIMITATIONS`; technical runs aggregated before inference |
| 4 | `PRJEB32839` | 150 | 75 | Full RNA-seq completed; `PASS_WITH_LIMITATIONS`; stage and sex contrasts |
| 5 | `E-MTAB-451` | 12 | 11 | Full RNA-seq completed; `PASS`; one technical-run pair aggregated before inference |
| 6 | `PRJEB3190` | 20 | 20 | Full RNA-seq completed; `PASS`; seven-point schistosomulum time course |
| 7 | `E-ERAD-478` | 120 | 60 | Full RNA-seq completed; `PASS`; two technical runs aggregated per biological sample |

Every study exposes the same compact QC interface at
`results/<study>/qc/qc_sample_summary.tsv` and `qc_summary.json`; detailed
MultiQC reports remain available alongside it where published. The fields,
classification rules, and generation procedure are defined in
[`docs/qc_results_contract.md`](docs/qc_results_contract.md).

`PRJNA602528` is the first completed real *S. mansoni* application and
operational calibration run. It is not a new benchmark of HelixForge. Its
public design has one library per time point, so the frozen project stops after
tximport and does not perform DESeq2 inference. The compact accepted products
and execution report are available under [`results/PRJNA602528/`](results/PRJNA602528/).

`PRJNA597909` is the first complete native RNA-seq application in this
project. All 20 paired-end samples passed QC, Salmon quantification and
tximport completed, and all four frozen DESeq2 contrasts were estimated. The
accepted tables, figures, portable terminal manifest, and report are available
under [`results/PRJNA597909/`](results/PRJNA597909/). Its
`PASS_WITH_LIMITATIONS` classification records terminal-runtime recovery and
the then-unresolved `-resume` incident; the scientific stages themselves
passed. Later investigation traced that incident to unsupported direct JAR
invocation rather than shared NFS, and the official launcher subsequently
passed complete cache-reuse and selective-invalidation validation.

`PRJEB14695` completed the native RNA-seq path for 138 paired-end technical
runs representing 23 biological samples. Run-to-sample mapping and aggregation
were enforced as scientific gates before tximport and DESeq2. The first attempt
exposed an overly strict validation of the run-level `technical_unit` field;
the reviewed HelixForge correction was then validated by a controlled resume
that reused all 815 eligible upstream tasks. The accepted matrices, six frozen
contrasts, figures, HTML reports, terminal manifest, and validation summary are
available under [`results/PRJEB14695/`](results/PRJEB14695/).

Earlier studies retain their recorded `v1.0.1` execution provenance.
`PRJEB14695` additionally pins the reviewed post-release contract fix at commit
`42864266892d1165477bb3b33c919e1fabb28ad1`; the current project default is
`v1.0.2`, as listed above, without rewriting prior runs.

`PRJEB32839` completed the `v1.0.2` native path for 150 technical runs and 75
biological samples. The exact frozen inventory, 12 DESeq2 contrasts, compact
QC/DE summaries, candidate-gene report record, terminal manifest, storage
calibration, and reviewed execution report are available under
[`results/PRJEB32839/`](results/PRJEB32839/). The newer pin applies to this
study and does not rewrite earlier execution provenance. Verified persistence
and the private audit archive preceded cleanup; 1.208 TB of study-specific
scratch data were removed without touching shared resources or other studies.

`E-MTAB-451` completed the native `v1.0.2` paired-end path for 12 runs
representing 11 biological samples. The two technical runs assigned to
`schistosomulum_3h_r1` were consolidated before quantification and inference.
All 136 processes completed, the three frozen `~ condition` contrasts were
estimated, and Salmon mapping rates ranged from 62.26% to 78.87% (median
75.18%). The accepted matrices, differential-expression results, QC summary,
terminal manifests, and sanitized execution reports are available under
[`results/E-MTAB-451/`](results/E-MTAB-451/). The complete candidate-gene
report passed but is retained outside Git under the current repository-size
policy. Persistence, checksums, private audit capture, and sensitive-path
scanning passed before the study-specific scratch data and FASTQs were removed.

`PRJEB3190` completed the native `v1.0.2` paired-end path for 20 biological
samples spanning seven schistosomulum time points. All 227 processes completed,
Salmon mapping rates ranged from 63.14% to 88.95% (median 82.95%), tximport
produced matrices for 9,914 genes, and all nine frozen `~ condition` contrasts
were estimated. ENA alias-derived sequencing labels remain exploratory because
they are not documented as technical batches and were not included in the
inferential model. The accepted matrices, differential-expression results, QC
summary, terminal manifests, and sanitized execution reports are available
under [`results/PRJEB3190/`](results/PRJEB3190/). The complete candidate-gene
report passed but remains outside Git under the repository-size policy.
Persistence, checksums, private audit capture, and sensitive-path scanning
passed before the study-specific workdir, FASTQs, and other reproducible
scratch intermediates were removed.

`E-ERAD-478` completed the native `v1.0.2` paired-end path for 120 runs
representing 60 biological samples. Two technical runs were aggregated per
sample before quantification and inference. All 950 processes completed,
Salmon mapping rates ranged from 75.03% to 91.84% (median 85.74%), tximport
produced matrices for 9,914 genes, and all 32 frozen contrasts were estimated
with the inferential design `~ batch + condition`. All 60 samples passed the
frozen QC review thresholds. The accepted matrices, differential-expression
results, QC summaries, terminal manifests, and sanitized execution reports are
available under [`results/E-ERAD-478/`](results/E-ERAD-478/). The complete
candidate-gene report passed but remains outside Git. Persistence, checksums,
private audit capture, and sensitive-path scanning passed before 621.87 GB of
study-specific scratch data were removed.

## Execution policy

Only one study may occupy heavy scratch storage at a time. Data acquisition is
outside the scientific Nextflow DAG, and all heavy work on an HPC system must
run inside Slurm allocations. Persistent results, terminal manifests,
provenance, and compact audit packages are retained; FASTQs, work directories,
indexes, and other reproducible intermediates are removed after acceptance.

See [server execution](docs/server_execution.md) and the
[server preparation report](docs/server_preparation_report.md) before staging
any data.

## RNA-seq atlas

The first integrated RNA-seq atlas combines seven accepted studies as a
descriptive, navigable report covering 219 biological samples, 470 technical
runs, 9,914 shared genes, and 66 original per-study DESeq2 contrasts. It
provides sample/QC inventory, global and within-study PCA, sample correlation,
selected-gene expression and effect maps, a searchable gene explorer, SVG
figures, machine-readable tables, and checksummed provenance. Canonical
`Smp_*` identifiers are accompanied by WBPS19 functional descriptions and
annotation provenance when available.

Open [`results/atlas/rnaseq/atlas.html`](results/atlas/rnaseq/atlas.html) or see
the [atlas contract and build instructions](analyses/atlas/README.md). The
global views are exploratory: the atlas does not pool studies for inference or
perform batch correction. Formal batch-effect assessment remains deferred.

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
