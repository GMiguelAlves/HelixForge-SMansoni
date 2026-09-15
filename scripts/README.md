# Operational scripts

- `download/`: reference and FASTQ acquisition through compute allocations.
- `validate/`: package, download, reference, runtime, and contract checks.
- `harmonize/`: conservative metadata indexing and controlled vocabularies.
- `summarize/`: accepted-run provenance collection.
- `run_study.sh`: HelixForge launcher for one registered study.

These helpers orchestrate project inputs; they do not replace HelixForge
scientific modules or perform nested Slurm submission from within Nextflow.

`validate/validate_reference.py` validates the materialized SM_V10/WBPS19
FASTA, annotations, tx2gene mapping and Salmon-index metadata. It reads large
reference files but does not modify them; execute it in a Slurm allocation on
the server. Site-specific `.sbatch` submission scripts are operational audit
material and are deliberately not versioned here.
