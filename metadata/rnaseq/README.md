# RNA-seq metadata

Each study directory contains `runs.tsv` (FASTQ acquisition and integrity),
`samples.tsv` (biological units), and `metadata.csv` (HelixForge input).
Where needed, `excluded_ena_runs.tsv` explicitly records public-deposit
exclusions.

`master_sample_manifest.tsv` is a conservative index of the initially
consolidated studies. `p1_pending_manifest.tsv` records a historical
preparation batch; it does not define a separate scientific project or
replace the individual study records. Files in `vocabularies/` preserve
curated values without assuming biological equivalence across studies.

Study contracts and current status are described in the
[study registry](../../docs/studies_and_preparation.md).
