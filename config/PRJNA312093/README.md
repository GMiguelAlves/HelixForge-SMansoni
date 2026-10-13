# PRJNA312093

## Scientific decision

Picard et al. studied sexual differentiation across development. Primary
source: https://doi.org/10.1371/journal.pntd.0004930.

The input package selected **20 RNA-seq runs** from *Schistosoma mansoni*,
representing **20 biological samples**. Layout: `SINGLE`. HelixForge mode:
`full`.

## Frozen limitation

The deposit combines 20 RNA-seq and 24 ChIP-seq libraries under one BioSample.
The design was reconstructed from library names. Fragment-length parameters
were not published and therefore required biological review.

The single-end fragment prior has mean `200` and standard deviation `80`.
These provisional protocol-based values were not reported by the paper and
required a sensitivity analysis before atlas inclusion.

Raw FASTQs are not versioned. `runs.tsv` freezes URLs, MD5 checksums, sizes,
and run-to-sample relationships. The immutable ENA snapshot is
`provenance/PRJNA312093/source_snapshots/ena_read_run.tsv`.
