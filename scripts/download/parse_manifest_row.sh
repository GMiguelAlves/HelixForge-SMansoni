#!/usr/bin/env bash

# Read the first 15 ENA manifest columns without collapsing empty TSV fields.
# The caller receives dataset, sample, run, layout, URLs, checksums and sizes.
download_parse_manifest_row() {
  local manifest=$1 task=$2
  local -a fields=()
  mapfile -d '' -t fields < <(
    awk -F '\t' -v row="$((task + 1))" '
      NR == 1 {
        if ($1 != "dataset" || $2 != "sample_id" || $3 != "run_accession" ||
            $9 != "library_layout" || $10 != "fastq_1_url" ||
            $11 != "fastq_2_url" || $12 != "fastq_1_md5" ||
            $13 != "fastq_2_md5" || $14 != "fastq_1_bytes" ||
            $15 != "fastq_2_bytes") exit 2
      }
      NR == row {
        if (NF < 15) exit 3
        for (i = 1; i <= 15; i++) printf "%s%c", $i, 0
        exit
      }
    ' "$manifest"
  )
  if [[ ${#fields[@]} -ne 15 ]]; then
    echo "Invalid ENA manifest header or no complete row for task $task: $manifest" >&2
    return 2
  fi
  dataset=${fields[0]}
  sample=${fields[1]}
  run=${fields[2]}
  biosample=${fields[3]}
  secondary=${fields[4]}
  experiment=${fields[5]}
  unit=${fields[6]}
  library=${fields[7]}
  layout=${fields[8]}
  url1=${fields[9]}
  url2=${fields[10]}
  md51=${fields[11]}
  md52=${fields[12]}
  bytes1=${fields[13]}
  bytes2=${fields[14]}
}
