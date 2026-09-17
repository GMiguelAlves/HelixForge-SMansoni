#!/usr/bin/env bash
set -euo pipefail
: "${HF_PACKAGE_ROOT:?source config/server.env}"
project=${1:?usage: submit_downloads.sh PROJECT}
settings="${HF_PACKAGE_ROOT}/config/${project}/user_settings.sh"
if [[ -s "$settings" ]]; then
  # shellcheck disable=SC1090
  source "$settings"
fi
manifest="${HF_PACKAGE_ROOT}/metadata/${project}/runs.tsv"
[[ -s "$manifest" ]]
count=$(($(wc -l < "$manifest") - 1))
mkdir -p "${HF_DATA_ROOT}/${project}/logs"
concurrency=${HF_DOWNLOAD_CONCURRENCY:-${HF_MAX_CONCURRENT_JOBS:-10}}
sbatch --job-name="dl_${project}"   --array="1-${count}%${concurrency}"   --cpus-per-task=1 --mem=2G --time=24:00:00   --output="${HF_DATA_ROOT}/${project}/logs/download_%A_%a.out"   --error="${HF_DATA_ROOT}/${project}/logs/download_%A_%a.err"   "${HF_PACKAGE_ROOT}/scripts/download/download_fastq_array.sh" "$manifest"
