# PRJNA704646 — descriptive RNA-seq execution review

## Scope and execution

[PRJNA704646](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA704646) contains
24 paired-end technical runs grouped into eight biological pools across
resistance status, adult/juvenile stage and sex. HelixForge used Salmon 1.10.3,
Nextflow 25.10.7 and the SM_V10/WBPS19 reference. All 48 input FASTQs passed
the checksum, gzip and pairing gate. The native QC → Salmon → Import path
completed; the trace has 171 successful tasks (169 `CACHED`, two `COMPLETED`).

The first Import attempt correctly rejected differing BioSample/library IDs
within pooled technical runs. The merged Import fix retains all such IDs when
each row declares its `technical_unit`, without relaxing biological metadata
conflict checks. The resumed run reused the eight Salmon quantifications; only
the corrected sample-table builder and tximport ran. No scientific parameter
or quantification result was changed by this metadata fix.

## Quality control

| Measure | Observed |
| --- | ---: |
| Biological pools | 8 |
| Technical runs | 24 |
| Pools classified `PASS` | 5/8 |
| Pools classified `REVIEW` | 3/8 |
| Minimum post-trim read retention | 96.52% |
| Salmon mapping range | 53.74–80.00% |
| MultiQC | `PASS` |
| Gene-level import | 9,914 genes × 8 pools |

The predefined review triggers are retention below 80% or Salmon mapping
below 60%. Only mapping triggered: `resistant_enriched_adult_male` (57.70%),
`sensitive_enriched_adult_female` (58.07%) and
`sensitive_enriched_adult_male` (53.74%). This is an operational review flag,
not evidence by itself of contamination or a reason to exclude samples. The
[QC summary](../qc/qc_summary.json) and
[sample table](../qc/qc_sample_summary.tsv) retain the exact values.

## Scientific boundary and provenance

The three runs per pool are technical observations of one biological unit;
they do not supply biological replication for a resistance, stage or sex
contrast within that pool. Therefore this execution intentionally stops at
gene-level counts, TPM and lengths. DESeq2 was not run. Any expression
comparison across the eight pools remains descriptive.

The native [Import manifest](../expression/import_manifest.json) records the
matrix checksums, all verified against this package. A terminal RNA-seq
Integration API manifest is not generated in `import` mode; the
[publication record](../manifests/publication_validation.json) says so
explicitly. A private archive preserves the full results, input metadata,
Nextflow trace, logs and checksums before study-scoped scratch cleanup. The
public MultiQC report removes 112 private work-directory paths but retains
the scientific plots and summary metrics.
