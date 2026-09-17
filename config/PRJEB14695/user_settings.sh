#!/usr/bin/env bash
: "${HF_PACKAGE_ROOT:?source config/server.env before running}"
: "${HF_DATA_ROOT:?set HF_DATA_ROOT}"
: "${HF_REFERENCE_ROOT:?set HF_REFERENCE_ROOT}"
: "${HF_RESULTS_ROOT:?set HF_RESULTS_ROOT}"

export PIPELINE_NAME="smansoni_prjeb14695_wbps19"
export ORGANISM_NAME="Schistosoma_mansoni"
export REFERENCE_ID="Schistosoma_mansoni_SM_V10_WBPS19"
export PIPELINE_PROJECTS="PRJEB14695"
export SCRATCH_ROOT="${HF_DATA_ROOT}"
export METADATA_FINAL="${HF_PACKAGE_ROOT}/metadata/PRJEB14695/metadata.csv"
export METADATA_FINAL_NEW="${METADATA_FINAL}"

export REF_GENOME_FA="${HF_REFERENCE_ROOT}/SM_V10_WBPS19/schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa"
export REF_TRANSCRIPTS_FA="${HF_REFERENCE_ROOT}/SM_V10_WBPS19/schistosoma_mansoni.PRJEA36577.WBPS19.mRNA_transcripts.fa"
export REF_GTF="${HF_REFERENCE_ROOT}/SM_V10_WBPS19/schistosoma_mansoni.PRJEA36577.WBPS19.canonical_geneset.gtf"
export REF_GFF3="${HF_REFERENCE_ROOT}/SM_V10_WBPS19/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3"
export SALMON_INDEX_DIR="${HF_REFERENCE_ROOT}/SM_V10_WBPS19/salmon_index_k31"
export QUANT_DIR="${HF_RESULTS_ROOT}/PRJEB14695/compatibility/quants"

export QUANT_METHOD="salmon"
export SALMON_KMER_SIZE=31
export TRIM_QUALITY=20
export TRIM_LENGTH=20
export RUN_BATCH_CORRECTION=0
export RUN_GENE_REPORT=0
export FASTQ_LAYOUT="paired"

# Project-specific operational ceiling agreed for the shared Slurm cluster.
export HF_MAX_CONCURRENT_JOBS=8
export HF_DOWNLOAD_CONCURRENCY=8
