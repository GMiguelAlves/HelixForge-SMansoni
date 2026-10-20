# Provenance

`source_package/` preserves the checksum manifest and review evidence used to
initialize this repository. Each RNA-seq study directory under `rnaseq/` is reserved for the exact
HelixForge pin, reference and annotation identities, command, environment,
Slurm metadata, checksums, terminal manifest, and execution summary produced by
its accepted run.

`reference/` records the verified SM_V10/WBPS19 checksums, compact validation
metrics and Salmon-index identity. `server_preparation_state.json` is the
sanitized, portable state of the initial server preflight; exact hostnames,
personal absolute paths and scheduler logs remain only in private audit
storage.

Older frozen records contain their original repository-relative paths. These
payloads and their checksums are preserved byte-for-byte; see
[`path_migration.md`](path_migration.md) to locate moved assay-specific files.
Future ChIP-seq execution records belong under `chipseq/`.
