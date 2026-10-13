# PRJEB14695

- Biological samples: **23**.
- Technical runs: **138**.
- Recommended mode: `full` (QC → Salmon → tximport → DESeq2).
- Consolidation: HelixForge receives one row per run, trims runs separately,
  and consolidates runs sharing a `sample_id` before Salmon.
- Historical preparation limit: no more than eight concurrent Slurm jobs for
  acquisition or workflow execution. Consult the current
  [server execution policy](../../docs/server_execution.md) before new runs.
- Statistical unit: DESeq2 receives 23 biological samples; the 138 technical
  runs are never treated as independent replicates.

## Groups

| Group | Biological samples |
|---|---:|
| `ovary_mixed_sex` | 3 |
| `ovary_single_sex` | 2 |
| `testis_mixed_sex` | 3 |
| `testis_single_sex` | 3 |
| `whole_female_mixed_sex` | 3 |
| `whole_female_single_sex` | 3 |
| `whole_male_mixed_sex` | 3 |
| `whole_male_single_sex` | 3 |

## Files and design provenance

- `metadata/rnaseq/PRJEB14695/runs.tsv`: ENA accessions, URLs, MD5, and bytes
  per run.
- `metadata/rnaseq/PRJEB14695/samples.tsv`: one row per biological sample.
- `metadata/rnaseq/PRJEB14695/metadata.csv`: run-level HelixForge input.
- `pipeline_config.sh` and `user_settings.sh`: portable adapters for the
  current input contract.
- `de_spec.json`: versioned design and contrasts, when estimable.

The launcher preserves the scientific fields of `de_spec.json` and creates
only an operational copy in the result directory. In that copy, relative
`target_dir` is resolved under the isolated project root. Checksums of the
source and operational copy are recorded in
`pipeline_info/runtime_specs/de_spec.provenance.json`.

The candidate-gene report is disabled for this study because no candidate
list was specified before analysis. This does not remove the HTML execution,
QC, and differential-expression reports required for review.
