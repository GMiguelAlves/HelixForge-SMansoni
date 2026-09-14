# Data dictionary

## Master manifest

`metadata/master_sample_manifest.tsv` contains one row per public sequencing
run. Values are copied conservatively from the curated study metadata; this
initialization does not collapse biological categories across studies.

| Field | Meaning |
|---|---|
| `study` | Public project accession. |
| `sample` | Stable project-local biological sample identifier. |
| `run` | ENA/SRA sequencing run accession. |
| `layout` | Library layout declared by ENA; current frozen studies are `PAIRED`. |
| `stage` | Curated life stage, preserved at study-specific resolution. |
| `sex` | Curated sex/composition label; mixed worms and worm pairs remain distinct. |
| `tissue` | Curated tissue or whole-organism sampling unit. |
| `condition` | Project-specific analysis group used by the frozen configuration. |
| `treatment` | Curated treatment/exposure label. |
| `time` | Time in hours or `not_applicable`, as supported by the study. |
| `replicate` | Biological replicate number within the project group. |
| `batch` | Supported batch label or `not_declared`; never inferred automatically. |
| `source` | Public metadata source, currently `ENA`. |

## Study files

- `metadata.csv`: run-level HelixForge input contract.
- `runs.tsv`: accession, library layout, download URL, MD5, byte size, and ENA
  provenance per technical run.
- `samples.tsv`: one row per biological sample and expected technical-run count.
- `excluded_ena_runs.tsv`: explicit exclusions where the live archive contains
  runs outside the frozen author-defined dataset.
- `config/<study>/de_spec.json`: reviewed DESeq2 design and contrasts.
- `config/PRJNA602528/de_spec.blocked.json`: explicit non-estimable DE state.

## Identifier conventions

`sample` and `run` identifiers are never inferred from filenames downstream.
Technical runs sharing a `sample` are trimmed individually and consolidated by
the HelixForge QC layer before Salmon. Reference identity is always
`Schistosoma_mansoni_SM_V10_WBPS19`; annotation identity is
`Schistosoma_mansoni_SM_V10_WBPS19_canonical_geneset`.

Controlled-vocabulary tables enumerate the exact source-preserving values used
at initialization. They are registries, not claims that similarly named stages,
tissues, sexes, or treatments are biologically interchangeable across studies.

