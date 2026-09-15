# Server execution policy

## Storage layout

Use environment variables rather than committed absolute paths:

- `HF_HELIXFORGE_ROOT`: checkout of HelixForge at the pinned release commit;
- `HF_PACKAGE_ROOT`: checkout of this repository;
- `HF_DATA_ROOT`: downloaded FASTQs on scratch;
- `HF_REFERENCE_ROOT`: shared SM_V10/WBPS19 reference and Salmon index;
- `HF_RESULTS_ROOT`: current study results;
- `HF_WORK_ROOT`: current study Nextflow work directory.

Copy `config/server.env.template` to the ignored `config/server.env`, set paths
for the target site, and source it. Never commit usernames, credentials,
private hostnames, home paths, scratch paths, SSH material, or tokens.

## Slurm rules

- The login/head node may host the Nextflow scheduler driver and light metadata
  validation only. Scientific processing and data transfer must run in Slurm
  allocations.
- Processes are submitted by Nextflow; scientific scripts must not submit
  nested jobs.
- Start with the site-safe queue limit in `config/slurm.config` and inspect
  cluster occupancy before increasing concurrency.
- Every compute node must see the same repository, reference, input, result,
  cache, and work paths required by the selected container profile.

## One-project-at-a-time policy

```text
ONE PROJECT AT A TIME
```

Only one study may hold heavy FASTQs, indexes under construction, work files,
and unpublished intermediates in scratch. The frozen operational order is:

1. `PRJNA602528` — Import-only operational calibration completed and cleaned;
2. `PRJNA597909`;
3. `PRJEB14695`;
4. `PRJEB32839`.

The order may change only with a documented operational reason. `PRJNA602528`
passed persistence and cleanup gates; starting another study still requires an
explicitly reviewed project-specific authorization.

## Checkpoints

1. Verify HelixForge tag/commit, Nextflow 25.10.7, Java 21, and the site profile.
2. Download and validate the reference in a compute allocation.
3. Download one study's FASTQs and verify ENA bytes, MD5, and gzip integrity.
4. Run the package and HelixForge preflight checks.
5. Execute the configured study mode and review trace, retries, QC, manifests,
   checksums, and scientific design status.
6. Copy compact accepted outputs and provenance to persistent storage.
7. Create and verify an audit package with a Portuguese README.
8. Remove only that study's verified reproducible intermediates from scratch.

`-resume` requires both the Nextflow task database and unchanged work outputs.
Its top-level persistence was environment-dependent during HelixForge v1
validation; confirm cache reuse on this site rather than assuming it.
During `PRJNA602528`, eligible QC tasks were not recovered from the shared NFS
task cache. A controlled post-QC re-entry consumed the already validated merged
FASTQs and completed the native Salmon and Import layers. This is an operational
runtime limitation and did not alter scientific parameters or results.

## Persistent records

Retain small QC summaries, DE tables, normalized counts/TPM when reasonably
sized, figures, HTML reports, terminal manifests with their portable artifact
bundles, performance summaries, commands, versions, checksums, Slurm metadata,
and execution summaries. Raw FASTQs, BAMs, BigWigs, indexes, work directories,
and large logs remain outside Git.
