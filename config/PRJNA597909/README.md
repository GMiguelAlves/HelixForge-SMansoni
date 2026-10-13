# PRJNA597909

- Biological samples: **20**.
- Technical runs: **20**.
- Recommended mode: `full` (QC → Salmon → tximport → DESeq2).
- Consolidation: HelixForge receives one row per run, trims runs separately,
  and consolidates runs sharing a `sample_id` before Salmon.

## Groups

| Group | Biological samples |
|---|---:|
| `adult_control` | 5 |
| `adult_pzq` | 5 |
| `juvenile_control` | 5 |
| `juvenile_pzq` | 5 |

## Files

- `metadata/rnaseq/PRJNA597909/runs.tsv`: ENA accessions, URLs, MD5, and bytes
  per run.
- `metadata/rnaseq/PRJNA597909/samples.tsv`: one row per biological sample.
- `metadata/rnaseq/PRJNA597909/metadata.csv`: run-level HelixForge input.
- `pipeline_config.sh` and `user_settings.sh`: portable adapters for the
  current input contract.
- `de_spec.json`: versioned design and contrasts, when estimable.
