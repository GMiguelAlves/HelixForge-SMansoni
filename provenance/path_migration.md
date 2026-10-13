# Repository path migration

The RNA-seq study files were grouped by assay without changing their contents.
For any RNA-seq study accession, the current repository path is obtained by
inserting `rnaseq/` after the first path component:

| Historical prefix | Current prefix |
|---|---|
| `config/<study>/` | `config/rnaseq/<study>/` |
| `results/<study>/` | `results/rnaseq/<study>/` |
| `provenance/<study>/` | `provenance/rnaseq/<study>/` |

This applies only to study directories. Shared files in `config/reference/`,
`provenance/reference/`, `provenance/source_package/`, and the roots of these
directories have not moved. Runtime paths under `HF_RESULTS_ROOT` on the Slurm
server are also unchanged.

Some frozen manifests and accepted result reports record the historical path.
Those fields are evidence of the original package and were deliberately not
rewritten or rehashed. Resolve them through this mapping when browsing the
current repository. The original tree remains available in Git history.
