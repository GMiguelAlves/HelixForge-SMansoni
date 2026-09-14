#!/usr/bin/env bash
set -euo pipefail
: "${HF_PACKAGE_ROOT:?source env/server.env}"
: "${HF_RESULTS_ROOT:?set HF_RESULTS_ROOT}"
dest="${HF_RESULTS_ROOT}/provenance_$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$dest"
cp "$HF_PACKAGE_ROOT/MANIFEST.sha256" "$dest/package_MANIFEST.sha256"
git -c safe.directory="$HF_HELIXFORGE_ROOT" -C "$HF_HELIXFORGE_ROOT" rev-parse HEAD > "$dest/helixforge_commit.txt"
java -version 2> "$dest/java_version.txt"
"${HF_NEXTFLOW_BIN:-$HF_HELIXFORGE_ROOT/nextflow}" -version > "$dest/nextflow_version.txt" 2>&1
cp "$HF_REFERENCE_ROOT/SM_V10_WBPS19/reference_checksums.sha256" "$dest/"
find "$HF_RESULTS_ROOT" -type f \( -name 'rnaseq_run_manifest.json' -o -name 'trace.tsv' -o -name 'execution_report.html' \) -print0 | sort -z | xargs -0 sha256sum > "$dest/result_artifacts.sha256"
echo "$dest"
