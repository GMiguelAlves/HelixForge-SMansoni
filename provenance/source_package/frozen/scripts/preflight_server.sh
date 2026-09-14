#!/usr/bin/env bash
set -euo pipefail
for key in HF_PACKAGE_ROOT HF_HELIXFORGE_ROOT HF_DATA_ROOT HF_REFERENCE_ROOT HF_RESULTS_ROOT HF_WORK_ROOT HF_HELIXFORGE_COMMIT; do [[ -n "${!key:-}" ]] || { echo "missing $key" >&2; exit 2; }; done
for cmd in git java sbatch squeue curl gzip md5sum sha256sum python3; do command -v "$cmd" >/dev/null || { echo "missing command: $cmd" >&2; exit 2; }; done
[[ $(git -c safe.directory="$HF_HELIXFORGE_ROOT" -C "$HF_HELIXFORGE_ROOT" rev-parse HEAD) == "$HF_HELIXFORGE_COMMIT" ]]
"${HF_NEXTFLOW_BIN:-$HF_HELIXFORGE_ROOT/nextflow}" -version
NEXTFLOW="${HF_NEXTFLOW_BIN:-$HF_HELIXFORGE_ROOT/nextflow}" "$HF_HELIXFORGE_ROOT/bin/helixforge-doctor"
test -s "$HF_REFERENCE_ROOT/SM_V10_WBPS19/REFERENCE_READY.json"
python3 "$HF_PACKAGE_ROOT/scripts/validate_package.py" "$HF_PACKAGE_ROOT"
echo 'SERVER_PREFLIGHT=PASS'
