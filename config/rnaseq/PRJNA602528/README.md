# PRJNA602528

- Biological samples: **10**.
- Technical runs: **10**.
- Recommended mode: `import` (QC → Salmon → tximport; no differential
  expression).
- Consolidation: HelixForge receives one row per run, trims runs separately,
  and consolidates runs sharing a `sample_id` before Salmon.

## Groups

| Group | Biological samples |
|---|---:|
| `time_0h` | 1 |
| `time_0p25h` | 1 |
| `time_12h` | 1 |
| `time_1h` | 1 |
| `time_24h` | 1 |
| `time_3h` | 1 |
| `time_48h` | 1 |
| `time_6h` | 1 |
| `time_96h` | 1 |
| `time_9h` | 1 |

The public deposit has one library per time point. A model with ten levels
would be saturated, with no residual replication. Therefore
`de_spec.blocked.json` is deliberately non-executable.

## Files

- `metadata/rnaseq/PRJNA602528/runs.tsv`: ENA accessions, URLs, MD5, and bytes
  per run.
- `metadata/rnaseq/PRJNA602528/samples.tsv`: one row per biological sample.
- `metadata/rnaseq/PRJNA602528/metadata.csv`: run-level HelixForge input.
- `pipeline_config.sh` and `user_settings.sh`: portable adapters for the
  current input contract.
- `de_spec.blocked.json`: explicit block on non-estimable inference.
