# PRJEB32839

- Biological samples: **75**.
- Technical runs: **150**.
- Recommended mode: `full` (QC → Salmon → tximport → DESeq2).
- Consolidation: HelixForge receives one row per run, trims runs separately,
  and consolidates runs sharing a `sample_id` before Salmon.
- Frozen runtime: HelixForge `v1.0.2`
  (`5d4b3e696319db5cd7633472504964f1dc7c0434`), without automatic
  updates from `master`.

## Groups

| Group | Biological samples |
|---|---:|
| `cercaria_female` | 5 |
| `cercaria_male` | 5 |
| `cercaria_mixed` | 5 |
| `egg_mixed` | 5 |
| `juvenile_26d_female` | 5 |
| `juvenile_26d_male` | 5 |
| `miracidium_mixed` | 5 |
| `schistosomulum_2d_female` | 5 |
| `schistosomulum_2d_male` | 5 |
| `schistosomulum_2d_mixed` | 5 |
| `sporocyst_1d_mixed` | 5 |
| `sporocyst_32d_female` | 5 |
| `sporocyst_32d_male` | 5 |
| `sporocyst_32d_mixed` | 5 |
| `sporocyst_5d_mixed` | 5 |

ENA currently lists 180 runs. This package freezes the 150 runs present in
the authors' official metadata; the remaining 30 are recorded in
`excluded_ena_runs.tsv`.

## Files

- `metadata/rnaseq/PRJEB32839/runs.tsv`: ENA accessions, URLs, MD5, and bytes
  per run.
- `metadata/rnaseq/PRJEB32839/samples.tsv`: one row per biological sample.
- `metadata/rnaseq/PRJEB32839/metadata.csv`: run-level HelixForge input.
- `pipeline_config.sh` and `user_settings.sh`: portable adapters for the
  current input contract.
- `de_spec.json`: versioned design and contrasts, when estimable.

## Storage gate

No FASTQ may be acquired before
`provenance/PRJEB32839/storage_preflight.json` records
`authorization: GO`. The calculation uses the frozen inventory of 150 runs,
the exact sum of 300 FASTQs, and a conservative 6× projection with a 20%
margin.
