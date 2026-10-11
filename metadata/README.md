# Metadata

Curated run- and sample-level records are versioned per study. The master
manifest is a conservative row-level index across the seven studies; controlled
vocabularies enumerate current values without performing biological
cross-study harmonization.

`p1_pending_manifest.tsv` is a separate preparation index for the seven P1
accessions not yet incorporated into the atlas. It does not modify the frozen
master manifest of completed studies. Each pending study contains `runs.tsv`
(download contract), `samples.tsv` (biological units) and `metadata.csv`
(HelixForge input). In particular, `PRJNA704646` consolidates three technical
libraries per biological RNA pool and `PRJNA312093` declares the native
single-end fragment prior required by HelixForge.
