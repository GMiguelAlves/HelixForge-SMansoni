# ChIP-seq metadata

This directory contains the public *Schistosoma mansoni* ChIP-seq input
contracts. It describes the files; the reviewed study-selection, control,
peak-type, execution, and atlas-admission decisions are maintained in the
[ChIP-seq core preparation record](../../docs/chipseq_core_preparation.md).
Presence in this directory does not by itself authorize analysis or atlas
inclusion.

## Shared registries

- `study_registry.tsv`: inclusion status, priority, executable scope, and
  raw-data volume by deposit.
- `run_inventory.tsv`: consolidated inventory of public FASTQs.
- `publication_registry.tsv`: deposit-to-publication links, including
  unresolved cases.
- `duplicate_registry.tsv` and `duplicate_map.tsv`: deposits and runs that
  must not be acquired or counted twice.

## Per-study files

- `runs.tsv`: ENA HTTPS URLs, MD5, byte sizes, and run counts for verifiable
  acquisition.
- `design.tsv`: curated stage, sex, treatment, biological unit, replicate,
  target, and control relationships, retaining public descriptions.
- `helixforge.narrow.tsv` and `helixforge.broad.tsv`: native HelixForge inputs
  separated by peak geometry.
- `helixforge.alignment.tsv`: QC/alignment-only input when an adequate public
  control is unavailable.

Raw FASTQs are not versioned. Local FASTQ names in `helixforge.*.tsv` match
the final path component of the corresponding URL in `runs.tsv`. The full ENA
query snapshots are retained under
`provenance/source_package/frozen/chipseq/<accession>/ena_read_run.tsv`.
