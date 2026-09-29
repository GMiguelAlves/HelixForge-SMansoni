#!/usr/bin/env bash
: "${HF_PACKAGE_ROOT:?source config/server.env}"
: "${HF_DATA_ROOT:?set HF_DATA_ROOT}"
: "${HF_REF_GENOME_FA:?set HF_REF_GENOME_FA}"
: "${HF_REF_TRANSCRIPTS_FA:?set HF_REF_TRANSCRIPTS_FA}"
: "${HF_REF_GFF3:?set HF_REF_GFF3}"
: "${HF_REF_GTF:?set HF_REF_GTF}"
: "${HF_SALMON_INDEX_DIR:?set HF_SALMON_INDEX_DIR}"
export PIPELINE_NAME="smansoni_prjeb3190"
export ORGANISM_NAME="Schistosoma mansoni"
export REFERENCE_ID="Schistosoma_mansoni_SM_V10_WBPS19"
export PIPELINE_PROJECTS="PRJEB3190"
export SCRATCH_ROOT="$HF_DATA_ROOT"
export METADATA_FINAL="$HF_PACKAGE_ROOT/metadata/PRJEB3190/metadata.csv"
export METADATA_FINAL_NEW="$METADATA_FINAL"
export REF_GENOME_FA="$HF_REF_GENOME_FA"
export REF_TRANSCRIPTS_FA="$HF_REF_TRANSCRIPTS_FA"
export REF_GFF3="$HF_REF_GFF3"
export REF_GTF="$HF_REF_GTF"
export SALMON_INDEX_DIR="$HF_SALMON_INDEX_DIR"
export QUANT_METHOD="salmon"
export FASTQ_LAYOUT="paired"
export PIPELINE_EXECUTOR="slurm"
export RUN_SALMON_INDEX=0
export RUN_STAR_INDEX=0
export RUN_STAR_GTF_INDEX=0
export RUN_BATCH_CORRECTION=0
export RUN_GENE_REPORT=0
export RNA_TOOLS_ENV="rna-tools"
export PYTHON_ENV="python-list"
export R_ANALYSIS_ENV="r-analysis"
