# Operational scripts

- `download/`: reference and FASTQ acquisition through compute allocations.
- `validate/`: package, download, reference, runtime, and contract checks.
- `harmonize/`: conservative metadata indexing and controlled vocabularies.
- `summarize/`: accepted-run provenance collection.
- `run_study.sh`: HelixForge launcher for one registered study.

These helpers orchestrate project inputs; they do not replace HelixForge
scientific modules or perform nested Slurm submission from within Nextflow.

