# RNA-seq QC results contract

Every registered RNA-seq study publishes the same two compact files under
`results/<study>/qc/`:

- `qc_sample_summary.tsv`: one row per biological sample;
- `qc_summary.json`: study-level aggregation and the rules used to classify
  samples.

MultiQC remains the detailed interactive report. These compact files are the
stable interface for documentation, review, cross-study summaries, and the UI.
They are derived from preserved FastQC and Salmon outputs; they do not replace
the original tool reports or alter any scientific result.

## Sample table

`qc_sample_summary.tsv` has the following columns in this order:

| Column | Meaning |
|---|---|
| `sample_id` | Biological sample identifier used by Import and DESeq2 |
| `raw_reads_both_mates` | Sum of FastQC `Total Sequences` over raw R1 and R2 files and all technical runs assigned to the sample |
| `trimmed_reads_both_mates` | Equivalent total after trimming |
| `trim_retention_percent` | `trimmed_reads_both_mates / raw_reads_both_mates * 100` |
| `salmon_mapping_percent` | Salmon `percent_mapped` for the biological sample |
| `classification` | `PASS`, `REVIEW`, or `FAIL` |
| `reason` | Semicolon-delimited machine-readable reasons; `none` for `PASS` |

The default review thresholds preserve the rules used in the accepted project
reports: trimming retention below 80% or Salmon mapping below 50%. A recorded
study-specific `REVIEW` or `FAIL` flag is preserved even when these two numeric
thresholds pass; it is never silently downgraded.

## Study summary

`qc_summary.json` uses schema version `1.0` and records the study identifier,
sample count, FastQC/MultiQC availability, total raw and retained reads,
retention distribution, Salmon mapping distribution, fragment totals, target
and software versions, classification counts, and the applied thresholds.

## Generation

Use the study-independent generator:

```bash
python3 scripts/summarize/build_qc_summary.py \
  --study-id STUDY \
  --results-root /path/to/preserved/results/STUDY \
  --output-dir results/STUDY/qc
```

For historical packages, `--audit-zip` reads the preserved archive directly.
For an already reviewed legacy sample table, `--sample-table` normalizes its
columns and preserves explicit review flags. New studies should use native
results directly and must not add another study-specific QC finalizer.
