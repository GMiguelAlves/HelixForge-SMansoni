#!/usr/bin/env bash
set -euo pipefail
: "${HF_DATA_ROOT:?source config/server.env}"
: "${HF_PACKAGE_ROOT:?source config/server.env}"
manifest=${1:?usage: download_fastq_array.sh metadata/rnaseq/PROJECT/runs.tsv}
task=${SLURM_ARRAY_TASK_ID:?submit as a Slurm array}
# shellcheck source=parse_manifest_row.sh
source "$HF_PACKAGE_ROOT/scripts/download/parse_manifest_row.sh"
download_parse_manifest_row "$manifest" "$task"
[[ "$layout" == "PAIRED" || "$layout" == "SINGLE" ]] || { echo "$run has unsupported layout: $layout" >&2; exit 2; }
[[ -n "$url1" && -n "$md51" && -n "$bytes1" ]] || { echo "$run has incomplete R1 metadata" >&2; exit 2; }
if [[ "$layout" == "PAIRED" ]]; then
  [[ -n "$url2" && -n "$md52" && -n "$bytes2" ]] || { echo "$run has incomplete R2 metadata" >&2; exit 2; }
else
  [[ -z "$url2" && -z "$md52" && -z "$bytes2" ]] || { echo "$run is single-end but R2 metadata is present" >&2; exit 2; }
fi
dest="${HF_DATA_ROOT}/${dataset}/fastq_ftp"
prov="${HF_DATA_ROOT}/${dataset}/download_provenance"
mkdir -p "$dest" "$prov"
download_one() {
  local url=$1 md5=$2 bytes=$3 output=$4
  if [[ -s "$output" ]] && [[ $(stat -c %s "$output") == "$bytes" ]] && echo "$md5  $output" | md5sum -c - >/dev/null 2>&1; then return 0; fi
  local part="${output}.part"
  curl --fail --location --retry 8 --retry-all-errors -C - -o "$part" "$url"
  [[ $(stat -c %s "$part") == "$bytes" ]]
  echo "$md5  $part" | md5sum -c -
  mv "$part" "$output"
}
download_one "$url1" "$md51" "$bytes1" "$dest/${run}_1.fastq.gz"
gzip -t "$dest/${run}_1.fastq.gz"
if [[ "$layout" == "PAIRED" ]]; then
  download_one "$url2" "$md52" "$bytes2" "$dest/${run}_2.fastq.gz"
  gzip -t "$dest/${run}_2.fastq.gz"
fi
printf '{"dataset":"%s","sample_id":"%s","run_accession":"%s","status":"DOWNLOAD_READY","checked_utc":"%s"}
' "$dataset" "$sample" "$run" "$(date -u +%FT%TZ)" > "$prov/${run}.json"
