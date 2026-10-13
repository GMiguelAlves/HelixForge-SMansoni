# PRJNA704646

## Scientific decision

Genetic basis of praziquantel response. Primary source:
https://doi.org/10.1126/scitranslmed.abj9114.

The input package selected **24 RNA-seq runs** from *Schistosoma mansoni*,
representing **8 biological samples**. Layout: `PAIRED`. HelixForge mode:
`import`.

## Frozen limitation

The three libraries for each combination came from the same RNA pool. They
are technical replicates and must be consolidated. After consolidation,
there is one biological unit per combination, so inferential differential
expression remains blocked.

Raw FASTQs are not versioned. `runs.tsv` freezes URLs, MD5 checksums, sizes,
and run-to-sample relationships. The immutable ENA snapshot is
`provenance/PRJNA704646/source_snapshots/ena_read_run.tsv`.
