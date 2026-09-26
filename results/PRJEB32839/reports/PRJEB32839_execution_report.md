# PRJEB32839 execution report

## Executive Summary

The frozen native RNA-seq analysis completed for 150 paired-end technical runs
representing 75 biological samples. Every scientific process completed: raw
and post-trim QC, technical-run merging, Salmon 1.10.3 quantification,
tximport, the frozen DESeq2 model, all 12 contrasts, a candidate-gene report,
and the terminal RNA-seq manifest. The main run executed 1,153 processes with
zero failures.

The final classification is `PASS_WITH_LIMITATIONS`. The limitations are
reported observations rather than execution failures: lower Salmon mapping in
all 15 32-day sporocyst samples, low trimming retention in six samples, 17
candidate IDs absent from the WBPS19 inputs, report generation through the
official re-entry mode, and the absence of continuous scratch peak sampling.
No sample, run, contrast, threshold, reference, or scientific parameter was
changed after results were observed.

## Frozen Dataset Selection

- Study: `PRJEB32839`
- Frozen biological samples: 75/75
- Frozen technical runs: 150/150
- Excluded public runs: 30, as preregistered
- Layout: paired-end
- FASTQ files: 300
- Duplicate assignments: 0
- Unmapped or orphan runs: 0
- Technical runs per biological sample: 2

The current public record is larger than the frozen project. Only the
preregistered 150-run set was admitted; no post-hoc expansion occurred.

## Storage Risk Assessment

| Quantity | Bytes |
|---|---:|
| Exact compressed FASTQ input | 258,443,852,490 |
| Central projection (5.007×) | 1,294,028,369,417 |
| Conservative projection (6.0×) | 1,550,663,114,940 |
| Extreme observed projection (6.606×) | 1,707,280,089,548 |
| Required capacity with 20% headroom | 1,860,795,737,928 |
| Scratch available at authorization | 2,269,043,752,960 |
| Capacity beyond required headroom | 408,248,015,032 |
| Warning threshold | 1,550,663,114,940 |
| Stop threshold | 2,155,591,565,312 |

No administrative per-project limit was provided; measured shared capacity
was therefore the effective preflight capacity, not an exclusive reservation.
The authorization gate was `GREEN/GO`. The observed final pre-cleanup footprint
was 1,207,940,984,325 bytes, below both thresholds.

## Dataset

The design contains 15 conditions, each with five biological replicates:
egg, miracidium, 1-day sporocyst, 5-day sporocyst, 32-day sporocyst, cercaria,
2-day schistosomulum, and 26-day juvenile stages, with mixed, female, and male
groups where defined by the frozen metadata.

## Biological Design

- Formula: `~ condition`
- Frozen factor: `condition`
- Conditions: 15
- Biological observations: 75
- Contrasts: 12
- Covariates: none
- Interactions: none
- Test: Wald
- Alpha: 0.05
- Absolute log2 fold-change threshold: 1.0
- Design matrix: full rank
- Requested contrasts: all estimable

Technical runs were trimmed independently and concatenated by mate before
Salmon. Only the resulting biological samples entered tximport and DESeq2, so
the analysis does not treat technical runs as biological replication.

## Reference

- Organism: *Schistosoma mansoni*
- Assembly: `SM_V10`
- Annotation: WormBase ParaSite `WBPS19`
- Reference ID: `Schistosoma_mansoni_SM_V10_WBPS19`
- Frozen reference checksums verified: 9
- Salmon index: prebuilt index reused; k=31; index version 5
- Index composite SHA-256: `88cf26aec437dc27959fa47afc8106f31943818e2260174796f0a0547223087b`
- Index reconstruction requested: no

## Acquisition

Acquisition used a checksum-aware Slurm array capped at eight concurrent jobs.
The initial pass completed 148 runs; one truncated transfer and one checksum
mismatch were isolated and retried. The two retries passed without changing
the frozen inventory. Acquisition and retries spanned 2026-09-24 18:08:44 to
2026-09-25 14:56:34 (20h47m50s).

## FASTQ Integrity

All 300 expected FASTQs were present. Official file sizes and MD5 values,
locally recorded SHA-256 values, gzip integrity, R1/R2 pairing, duplicate-name
checks, and orphan checks passed. Observed compressed bytes exactly matched
the preregistered 258,443,852,490 bytes.

## Native Execution

- HelixForge release: `v1.0.2`
- HelixForge commit: `5d4b3e696319db5cd7633472504964f1dc7c0434`
- Nextflow: `25.10.7` build 12755
- Java: Temurin 21.0.12+8
- Python: 3.10.18
- jsonschema: 4.25.1
- Salmon: 1.10.3
- Run name: `prjeb32839-full-v102-20260925`
- Main workflow wall time: 8h52m29s
- CPU time reported by Nextflow: 332.1 hours
- Processes: 1,153 completed; 0 failed
- Scheduler concurrency cap: 8
- Resume required: no

The certified Python runtime was first in `PATH`, `jsonschema` imported before
science, and compute jobs did not depend on Git. A failed administrative
validation attempt used an obsolete Python path and exited before reading any
scientific artifact; the certified runtime was then used. This did not alter
the workflow or its products.

## QC

- FastQC HTML reports: 750
- MultiQC: PASS
- Raw reads across both mates: 8,816,561,760
- Retained reads across both mates: 8,144,155,728
- Mean trimming retention: 92.28%
- Median trimming retention: 95.08%
- Trimming retention range: 41.52%–98.80%
- Salmon samples: 75
- Processed fragments: 4,072,077,864
- Mean Salmon mapping: 69.89%
- Median Salmon mapping: 80.36%
- Salmon mapping range: 17.25%–86.96%
- QC classifications: 54 PASS, 21 REVIEW, 0 FAIL

Six samples were marked `REVIEW` for trimming retention below 80%. All 15
32-day sporocyst samples were marked `REVIEW` for Salmon mapping below 50%, a
coherent stage-specific pattern rather than isolated random failures. Samples
were not removed automatically. The full per-sample classification and reason
are recorded in `qc_sample_summary.tsv`.

The frozen outputs do not contain a whole-transcriptome PCA, sample-correlation
matrix, or distance matrix. Candidate-gene PCA/MDS plots exist, but they are
not substituted for whole-transcriptome outlier diagnostics. This is recorded
as a reporting limitation and did not trigger post-hoc filtering.

## Expression

- Gene count matrix: 9,914 genes × 75 samples
- TPM matrix: 9,914 genes × 75 samples
- Effective-length matrix: 9,914 genes × 75 samples
- Normalized counts: preserved
- Salmon transcript targets per sample: 10,913
- `countsFromAbundance`: `lengthScaledTPM`
- `ignoreTxVersion=TRUE`
- `ignoreAfterBar=TRUE`

## Differential Expression

Each frozen contrast tested 9,792 genes. Significant genes satisfy
`padj < 0.05` and `abs(log2FoldChange) >= 1`. Positive fold change denotes
higher expression in the numerator.

| Contrast | Significant | Up | Down |
|---|---:|---:|---:|
| Miracidium mixed vs egg mixed | 3,394 | 2,061 | 1,333 |
| Sporocyst 1d mixed vs miracidium mixed | 3,585 | 868 | 2,717 |
| Sporocyst 5d mixed vs sporocyst 1d mixed | 1,681 | 739 | 942 |
| Sporocyst 32d mixed vs sporocyst 5d mixed | 2,942 | 2,048 | 894 |
| Cercaria mixed vs sporocyst 32d mixed | 4,423 | 2,000 | 2,423 |
| Schistosomulum 2d mixed vs cercaria mixed | 3,429 | 1,726 | 1,703 |
| Juvenile 26d female vs schistosomulum 2d female | 2,622 | 1,619 | 1,003 |
| Juvenile 26d male vs schistosomulum 2d male | 2,635 | 1,637 | 998 |
| Sporocyst 32d female vs male | 10 | 6 | 4 |
| Cercaria female vs male | 519 | 280 | 239 |
| Schistosomulum 2d female vs male | 680 | 366 | 314 |
| Juvenile 26d female vs male | 662 | 320 | 342 |

Independent filtering and complete effect distributions remain preserved in
the DESeq2 outputs. These results were not used to alter the frozen design.

## Candidate-gene Report

The official `report_reentry` mode consumed the preserved Import and DE
artifacts and ran only report context plus report generation. It completed in
4m49s without recomputing scientific processes. The report contains 65
preregistered queries in four groups, 48 matched genes, 17 explicit unmatched
identifiers, five audit tables, one HTML report, and 348 figures. The full
linked report tree is retained in persistent storage; the public bundle keeps
compact tables and overview figures.

## Performance

The execution demonstrates operational feasibility in Slurm, not hardware
scaling. The concurrency cap was eight jobs. FastQC had the largest observed
task RSS (approximately 1.50 GiB); Salmon reached approximately 1.10 GiB and
the largest DESeq2 contrast approximately 1.30 GiB. Per-process counts, summed
runtime, maximum runtime, and maximum RSS are recorded in `performance.tsv`.

## Scratch Calibration

| Study | Runs | Samples | FASTQ bytes | Observed footprint | Multiplier | Measurement |
|---|---:|---:|---:|---:|---:|---|
| PRJNA602528 | 10 | 10 | 17,455,857,588 | 97,859,705,391 | 6.606× | conservative observed footprint |
| PRJNA597909 | 20 | 20 | 33,199,877,299 | 166,370,873,344 | 5.007× | clean observed footprint |
| PRJEB14695 | 138 | 23 | 89,131,626,404 | 417,136,095,414 | 4.680× | final pre-cleanup footprint |
| PRJEB32839 | 150 | 75 | 258,443,852,490 | 1,207,940,984,325 | 4.674× | final pre-cleanup footprint |

PRJEB32839 was measured after download and immediately before cleanup, but not
at every requested logical checkpoint. Therefore 1.208 TB and 4.674× are
observed final pre-cleanup values, not claims of an exact continuous peak.
Neither the warning nor stop threshold was observed as crossed.

## Terminal Manifest

The terminal RNA-seq manifest passed JSON Schema, semantic, and tracked-input
filesystem validation. It declares 75 transcript-abundance artifacts, one
gene-count matrix, one gene-abundance matrix, one normalized-count matrix, 12
differential-expression artifacts, and one aggregate DE summary.

## Provenance

Runtime versions, reference and index checksums, frozen run inventory, official
and local FASTQ checksums, design matrix, contrast orientation, per-process
performance, per-sample QC, DE summaries, report manifests, and terminal
manifest validation are retained. Public records are path-sanitized; complete
private operational evidence is stored separately in the verified audit
archive.

Audit archive: `helixforge-smansoni-PRJEB32839-2026-09-26.zip`  
Audit SHA-256: recorded in the adjacent `.zip.sha256` file and in the final
cleanup provenance record.

## Limitations

- Fifteen 32-day sporocyst samples mapped below 50%; the coherent stage-level
  pattern requires biological review and was not used for automatic exclusion.
- Six samples retained below 80% of raw reads after trimming; they remain
  explicit `REVIEW` observations.
- Whole-transcriptome PCA, sample-correlation, and distance plots were not
  produced by the frozen path; candidate-gene ordinations are not treated as a
  substitute.
- Seventeen preregistered candidate identifiers were absent from WBPS19 and
  remain explicit unmatched records.
- Reporting was disabled in the full invocation and completed through the
  supported report re-entry mode; no scientific process was recomputed.
- Storage was not sampled continuously, so the exact peak and exact peak
  multiplier remain unavailable.

## Cleanup

Cleanup ran only after the persistent result copy, audit ZIP, checksums,
report, and public sanitization had passed. It recovered 1,207,941,033,218
bytes from the study-specific scratch root and left zero bytes of study data
or heavy residuals there. Raw and trimmed FASTQs, the Nextflow work directory,
temporary outputs, and study-local runtime state were removed. Shared
SM_V10/WBPS19 resources, the validated Salmon index, Nextflow runtime/cache,
persistent results, and all other studies were outside the cleanup scope and
were not touched.

## Final Classification

```text
FROZEN_RUN_INVENTORY = PASS
STORAGE_PREFLIGHT = PASS
STORAGE_AUTHORIZATION = PASS
DATA_ACQUISITION = PASS
FASTQ_INTEGRITY = PASS
POST_DOWNLOAD_STORAGE_GATE = PASS
REFERENCE_INTEGRITY = PASS
SALMON_INDEX_REUSE = PASS
PYTHON_RUNTIME_PREFLIGHT = PASS
METADATA_VALIDATION = PASS
STATISTICAL_DESIGN = PASS
REPORT_CONTRACT_PREFLIGHT = PASS
TECHNICAL_EXECUTION = PASS
QC = PASS_WITH_LIMITATIONS
EXPRESSION_IMPORT = PASS
DIFFERENTIAL_EXPRESSION = PASS
TERMINAL_MANIFEST = PASS
PROVENANCE = PASS
SCRATCH_CALIBRATION = PASS_WITH_LIMITATIONS
PERSISTENT_COPY = PASS
AUDIT_PACKAGE = PASS
CLEANUP = PASS

RESUME_REQUIRED = NO
RESUME_BEHAVIOR = NOT_APPLICABLE
REENTRY_REQUIRED = YES

PRJEB32839_RNASEQ_ANALYSIS = PASS_WITH_LIMITATIONS
```

`READY_FOR_PRJEB32839_REVIEW`
