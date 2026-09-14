#!/usr/bin/env bash
set -euo pipefail
: "${HF_DATA_ROOT:?source config/server.env}"
manifest=${1:?usage: download_fastq_array.sh metadata/PROJECT/runs.tsv}
task=${SLURM_ARRAY_TASK_ID:?submit as a Slurm array}
line=$(awk -F '	' -v n="$((task + 1))" 'NR==n {print; exit}' "$manifest")
[[ -n "$line" ]] || { echo "No manifest row for task $task" >&2; exit 2; }
IFS=$'	' read -r dataset sample run biosample secondary experiment unit library layout url1 url2 md51 md52 bytes1 bytes2 rest <<< "$line"
[[ "$layout" == "PAIRED" ]] || { echo "$run is not paired" >&2; exit 2; }
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
download_one "$url2" "$md52" "$bytes2" "$dest/${run}_2.fastq.gz"
gzip -t "$dest/${run}_1.fastq.gz"
gzip -t "$dest/${run}_2.fastq.gz"
printf '{"dataset":"%s","sample_id":"%s","run_accession":"%s","status":"DOWNLOAD_READY","checked_utc":"%s"}
' "$dataset" "$sample" "$run" "$(date -u +%FT%TZ)" > "$prov/${run}.json"
