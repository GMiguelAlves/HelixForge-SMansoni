# E-ERAD-478 results

The complete native paired-end RNA-seq execution is classified `PASS`. The
workflow completed 950/950 processes without failures and reused the certified
WBPS19/SM_V10 Salmon index.

- 120 paired-end ENA runs were aggregated into 60 biological samples, with two
  technical runs per sample.
- Salmon produced 60 sample quantifications; all 60 samples passed the frozen
  QC review thresholds.
- tximport produced gene-level count, TPM and effective-length matrices with
  9,914 genes and 60 sample columns.
- DESeq2 used the frozen `~ batch + condition` design and completed all 32
  prespecified contrasts across sex, pairing state and collection day.
- MultiQC, the candidate-gene report, Nextflow reports and the terminal RNA-seq
  manifest completed successfully. The full candidate-gene report and the
  re-entry package remain in verified private project storage.

Directory guide:

- `qc/`: standardized sample-level QC tables.
- `expression/`: gene-level matrices, import statistics and R provenance.
- `differential_expression/`: normalized counts, contrast tables and figures.
- `reports/`: MultiQC, Nextflow reports, validation and performance evidence.
- `manifests/`: terminal, reference and quantification manifests.

The public package intentionally excludes complete Salmon directories,
duplicated Integration API payloads, the RDS object and the complete
candidate-gene report. Those artifacts are preserved in the verified private
result and audit packages. This keeps the Git tree reviewable without changing
the scientific result or the re-entry evidence.

The authoritative terminal manifest is
[`manifests/rnaseq_run_manifest.json`](manifests/rnaseq_run_manifest.json).
The compact Git package is a publication view; its terminal manifest records
provenance, while full re-entry is validated against the preserved audit
package.

`E_ERAD_478_RNASEQ_ANALYSIS = PASS`
