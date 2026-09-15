# Server preparation report

## Scope

This report records the infrastructure and reference preparation completed on
2026-09-14 before the first `PRJNA602528` run. No study FASTQ was downloaded
and no scientific HelixForge workflow was executed.

## Environment

| Check | Frozen value | Status |
|---|---|---|
| HelixForge | `v1.0.0` | PASS |
| HelixForge commit | `14dc5a75d6c63f20d10c135f0c75138ea76dcc12` | PASS |
| Nextflow | `25.10.7` | PASS |
| Java | `21.0.12` | PASS |
| Salmon | `1.10.3` | PASS |
| R / rtracklayer | `4.3.3` / `1.62.0` | PASS |

The HelixForge checkout is detached at the release commit so future analysis
cannot silently follow `master`. The project checkout was clean before the
preparation changes. Public records use symbolic paths; personal server paths
and the head-node name remain only in the private audit record.

## Slurm

The institutional Slurm cluster exposed four compute nodes through the
`general` partition. `sbatch`, `squeue`, `sacct`, and `sinfo` were available.
A minimal one-CPU sanity allocation confirmed the pinned Java, Nextflow and
Salmon runtimes, shared mounts and compute-node memory. The job ran for three
seconds and started no scientific process.

The project profile limits Nextflow to five queued jobs and ten submissions per
second. Nextflow remains the only scheduler of scientific processes; nested
submission is prohibited.

## Storage

| Measurement before `PRJNA602528` | Value |
|---|---:|
| Shared home filesystem available | 702,530,191,360 bytes |
| Shared scratch filesystem available | 6,299,499,495,424 bytes |
| User home content | 35,087,294,172 bytes |
| User scratch content before project | 3 bytes |
| Prepared reference directory | 260,553,475 bytes |
| Reused Salmon index | 158,307,779 bytes |
| Project checkout | 1,085,747 bytes |

The operational warning threshold is 100 GiB and the stop threshold is 150
GiB. `ONE_PROJECT_AT_A_TIME = TRUE`. Software, Nextflow state, compact logs and
audit material remain in persistent storage; study FASTQs and the independent
work directory will use scratch. The reference and index are persistent shared
resources, not study intermediates.

## SM_V10 and WBPS19

The authority remains exactly:

- organism: *Schistosoma mansoni*;
- assembly: `SM_V10`;
- source/release: WormBase ParaSite `WBPS19`;
- reference ID: `Schistosoma_mansoni_SM_V10_WBPS19`.

The four official compressed artifacts were retrieved from the frozen URLs,
passed gzip integrity checks and received SHA-256 values. Existing V10 genome,
GFF3 and Salmon index resources were reused only after identity and format
verification. The official transcript FASTA and canonical GTF were materialized
from the verified WBPS19 files. This reuse avoids rebuilding an identical index
without changing the scientific reference.

Full checksums are in
[`provenance/reference/reference_checksums.sha256`](../provenance/reference/reference_checksums.sha256).

## Reference and annotation validation

| Artifact/check | Result |
|---|---:|
| Genome sequences / bases | 10 / 391,422,299 |
| Transcript sequences / bases | 10,960 / 39,461,259 |
| GFF3 records / genes / transcripts | 1,500,427 / 9,920 / 10,935 |
| Canonical GTF records / genes / transcripts | 241,537 / 9,920 / 10,960 |
| Duplicate FASTA identifiers | 0 |
| Invalid FASTA characters | 0 |
| Annotation contigs absent from genome | 0 |
| Reference-annotation compatibility | PASS |

WBPS19 identifiers are preserved as supplied. No suffix stripping, arbitrary
renaming or ID normalization was introduced during preparation.

## tx2gene

The mapping was produced with the HelixForge `v1.0.0` production policy from
the canonical GTF:

| Metric | Value |
|---|---:|
| Rows / transcripts | 10,960 |
| Genes | 9,920 |
| Transcriptome IDs mapped | 10,960 |
| Unmapped transcriptome IDs | 0 |
| Duplicate pairs | 0 |
| Transcripts assigned to multiple genes | 0 |

## Salmon index

The persistent index was originally built with Salmon `1.10.3` and `k=31` and
was accepted after validating its metadata, inventory, reference identity and
compatibility. It contains 15 files, occupies 158,307,779 bytes and has
composite SHA-256
`88cf26aec437dc27959fa47afc8106f31943818e2260174796f0a0547223087b`.
The reproducible command template and reuse decision are recorded in the
[index manifest](../provenance/reference/salmon_index_manifest.json).

## PRJNA602528 preflight

| Contract | Result |
|---|---|
| Runs | 10 |
| Samples | 10 |
| Layout | paired-end |
| R1/R2 completeness | PASS |
| Duplicate accessions | 0 |
| ENA FASTQ URL syntax | PASS |
| ENA HEAD requests | 20/20 HTTP 200 |
| Project mode | `IMPORT_ONLY` |
| Differential expression | `BLOCKED` |
| FASTQs downloaded | no |
| Scientific workflow executed | no |

Differential expression remains blocked because the publicly resolved
replication structure is insufficient for the frozen differential-expression
design. This is a property of the study design, not a HelixForge limitation.

## Permissions and paths

Compute-node probes passed for reference reads and writes to the future raw,
work, result, log, state and persistent audit destinations. No permission probe
file remained afterward. The operational config is private and ignored by Git;
tracked project files contain no personal absolute server path, Windows path,
credential, token or SSH material.

## Validation and CI

The reusable validator, compact contract tests, JSON/TSV checks, shell syntax
checks, path/secrets scan and Nextflow configuration parsing passed. The
generic discovery command executed all 10 tests successfully, and the GitHub
Actions `contracts` job passed for the preparation pull request.

## Limitations

- The site does not expose a user quota command, so the documented 150 GiB
  home quota is an administrative limit rather than a command-derived value.
- Git is available on the head node but not on compute nodes. The pinned commit
  is verified before submission; compute jobs validate shared artifacts and
  checksums without adding a dedicated Git environment.
- Available filesystem figures are shared NFS values, not a private allocation.
- This is an operational preparation record, not a performance benchmark.
- Nextflow `-resume` will be evaluated during a later controlled execution; no
  cache-reuse claim is made here.

## Readiness decision

All frozen scientific and operational preconditions are satisfied. The server
is ready for the separately authorized acquisition and import-only execution
of `PRJNA602528`. This report does not itself authorize or start that work.

```text
HELIXFORGE_ENVIRONMENT = PASS
SLURM_ENVIRONMENT = PASS
STORAGE_PREFLIGHT = PASS
SM_V10_REFERENCE = PASS
WBPS19_ANNOTATION = PASS
REFERENCE_ANNOTATION_COMPATIBILITY = PASS
TX2GENE_RESOURCE = PASS
SALMON_INDEX = PASS
PRJNA602528_METADATA = PASS
PRJNA602528_CONFIG = PASS
PRJNA602528_IMPORT_ONLY_CONTRACT = PASS
TEST_SUITE = PASS
CI = PASS

SERVER_PREPARATION = PASS
```
