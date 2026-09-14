#!/usr/bin/env bash
set -euo pipefail
: "${HF_PACKAGE_ROOT:?source config/server.env}"
: "${HF_HELIXFORGE_ROOT:?set HF_HELIXFORGE_ROOT}"
: "${HF_RESULTS_ROOT:?set HF_RESULTS_ROOT}"
: "${HF_WORK_ROOT:?set HF_WORK_ROOT}"
project=${1:?usage: run_study.sh PROJECT}
case "$project" in
  PRJEB32839|PRJEB14695|PRJNA597909) mode=full ;;
  PRJNA602528) mode=import ;;
  *) echo "Unsupported project: $project" >&2; exit 2 ;;
esac
current=$(git -c safe.directory="$HF_HELIXFORGE_ROOT" -C "$HF_HELIXFORGE_ROOT" rev-parse HEAD)
[[ "$current" == "$HF_HELIXFORGE_COMMIT" ]] || { echo "HelixForge commit mismatch: $current" >&2; exit 2; }
nextflow_bin=${HF_NEXTFLOW_BIN:-$HF_HELIXFORGE_ROOT/nextflow}
command -v "$nextflow_bin" >/dev/null 2>&1 || [[ -x "$nextflow_bin" ]]
config="$HF_PACKAGE_ROOT/config/$project/pipeline_config.sh"
outdir="$HF_RESULTS_ROOT/$project"
args=(run "$HF_HELIXFORGE_ROOT" -profile "${HF_PROFILES:-apptainer,slurm}" -c "$HF_PACKAGE_ROOT/config/slurm.config" -w "$HF_WORK_ROOT" -resume --workflow rnaseq --rnaseq_run_mode "$mode" --rnaseq_analysis_mode quantification --rnaseq_config "$config" --rnaseq_import_policy production_v1 --rnaseq_library_protocol full_length --rnaseq_counts_from_abundance lengthScaledTPM --salmon_validate_mappings true --rnaseq_report_enabled false --outdir "$outdir" -with-report "$outdir/execution_report.html" -with-trace "$outdir/trace.tsv" -with-timeline "$outdir/timeline.html" -with-dag "$outdir/dag.html")
if [[ "$mode" == full ]]; then args+=(--rnaseq_de_spec "$HF_PACKAGE_ROOT/config/$project/de_spec.json"); fi
"$nextflow_bin" "${args[@]}"
