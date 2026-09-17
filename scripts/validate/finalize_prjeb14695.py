#!/usr/bin/env python3
"""Validate and summarize the completed PRJEB14695 HelixForge run."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


EXPECTED_SAMPLES = 23
EXPECTED_RUNS = 138
EXPECTED_FASTQS = 276
EXPECTED_GENES = 9914
EXPECTED_TESTED_GENES = 9517
EXPECTED_CONTRASTS = 6
EXPECTED_FASTQC_HTML = 598


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def matrix_shape(path: Path) -> tuple[int, int, list[str], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        genes = [row[0] for row in reader]
    return len(genes), len(header) - 1, header[1:], genes


def seconds(value: str) -> float:
    if not value or value == "-":
        return 0.0
    units = {"ms": 0.001, "s": 1.0, "m": 60.0, "h": 3600.0, "d": 86400.0}
    return sum(float(number) * units[unit] for number, unit in re.findall(r"([0-9.]+)\s*(ms|s|m|h|d)", value))


def memory_bytes(value: str) -> int:
    if not value or value == "-":
        return 0
    match = re.fullmatch(r"\s*([0-9.]+)\s*([KMGTP]?B?|bytes?)\s*", value, re.IGNORECASE)
    if not match:
        return 0
    number = float(match.group(1))
    unit = match.group(2).upper().rstrip("B")
    factors = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4, "P": 1024**5}
    return int(number * factors.get(unit, 1))


def require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--fastq-integrity", required=True, type=Path)
    parser.add_argument("--sample-accounting", required=True, type=Path)
    parser.add_argument("--helixforge-commit", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    results = args.results
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    errors: list[str] = []

    integrity = json.loads(args.fastq_integrity.read_text(encoding="utf-8"))
    accounting = json.loads(args.sample_accounting.read_text(encoding="utf-8"))
    require(integrity.get("status") == "PASS", "FASTQ integrity did not pass", errors)
    require(integrity.get("runs") == EXPECTED_RUNS, "unexpected run count", errors)
    require(integrity.get("observed_fastq_files") == EXPECTED_FASTQS, "unexpected FASTQ count", errors)
    require(accounting.get("status") == "PASS", "sample accounting did not pass", errors)
    require(accounting.get("biological_samples") == EXPECTED_SAMPLES, "unexpected biological sample count", errors)
    require(accounting.get("deseq2_observations") == EXPECTED_SAMPLES, "DESeq2 observations are not biological samples", errors)

    import_root = results / "pipeline_info/native_import"
    matrix_root = import_root / "tximport"
    matrices = {}
    matrix_samples: list[str] | None = None
    matrix_genes: list[str] | None = None
    for name, filename in (("counts", "counts_matrix.tsv"), ("tpm", "tpm_matrix.tsv"), ("length", "length_matrix.tsv")):
        rows, columns, samples, genes = matrix_shape(matrix_root / filename)
        matrices[name] = {"genes": rows, "samples": columns}
        require(rows == EXPECTED_GENES, f"{name} matrix has {rows} genes", errors)
        require(columns == EXPECTED_SAMPLES, f"{name} matrix has {columns} samples", errors)
        if matrix_samples is None:
            matrix_samples, matrix_genes = samples, genes
        else:
            require(samples == matrix_samples, f"{name} sample order differs", errors)
            require(genes == matrix_genes, f"{name} gene order differs", errors)

    sample_rows = read_tsv(import_root / "sample_tables/import_samples.tsv")
    sample_ids = [row["sample_id"] for row in sample_rows]
    require(len(sample_rows) == EXPECTED_SAMPLES, "import sample table does not contain 23 rows", errors)
    require(len(sample_ids) == len(set(sample_ids)), "duplicated biological sample in import table", errors)

    de_root = results / "pipeline_info/native_de"
    de_preflight = json.loads((de_root / "preflight/preflight_report.json").read_text(encoding="utf-8"))
    require(de_preflight.get("status") == "valid", "DE preflight is not valid", errors)
    require(de_preflight.get("samples") == EXPECTED_SAMPLES, "DE preflight sample count differs", errors)
    require(de_preflight.get("genes") == EXPECTED_GENES, "DE preflight gene count differs", errors)
    require(de_preflight.get("contrasts") == EXPECTED_CONTRASTS, "DE preflight contrast count differs", errors)

    contrast_rows = read_tsv(de_root / "aggregate/deg_summary.tsv")
    require(len(contrast_rows) == EXPECTED_CONTRASTS, "unexpected number of completed contrasts", errors)
    require(len({row["contrast"] for row in contrast_rows}) == EXPECTED_CONTRASTS, "duplicated contrast result", errors)
    for row in contrast_rows:
        require(row["status"] == "ok", f"contrast {row['contrast']} did not pass", errors)
        require(int(row["n_samples"]) == EXPECTED_SAMPLES, f"contrast {row['contrast']} has wrong sample count", errors)
        require(int(row["n_genes"]) == EXPECTED_TESTED_GENES, f"contrast {row['contrast']} has wrong tested-gene count", errors)

    quant_root = results / "pipeline_info/native_quantification/salmon_quant"
    meta_files = sorted(quant_root.glob("*.quantification.quantification_statistics/meta_info.json"))
    require(len(meta_files) == EXPECTED_SAMPLES, "missing Salmon meta_info records", errors)
    mapping = []
    targets = set()
    salmon_versions = set()
    for path in meta_files:
        document = json.loads(path.read_text(encoding="utf-8"))
        mapping.append(float(document["percent_mapped"]))
        targets.add(int(document["num_valid_targets"]))
        salmon_versions.add(str(document["salmon_version"]))
    require(mapping and all(math.isfinite(value) and 0 <= value <= 100 for value in mapping), "invalid Salmon mapping percentages", errors)
    require(targets == {10913}, f"unexpected Salmon target counts: {sorted(targets)}", errors)
    require(salmon_versions == {"1.10.3"}, f"unexpected Salmon versions: {sorted(salmon_versions)}", errors)

    fastqc_count = len(list((results / "pipeline_info/native_qc/fastqc").glob("*.html")))
    multiqc = results / "pipeline_info/native_qc/multiqc/PRJEB14695_multiqc_030.html"
    require(fastqc_count == EXPECTED_FASTQC_HTML, f"unexpected FastQC HTML count: {fastqc_count}", errors)
    require(multiqc.is_file() and multiqc.stat().st_size > 0, "MultiQC HTML is missing", errors)

    validation = json.loads((results / "pipeline_info/integration_api/run_manifest.validation.json").read_text(encoding="utf-8"))
    terminal_manifest = json.loads((results / "pipeline_info/integration_api/rnaseq_run_manifest.json").read_text(encoding="utf-8"))
    require(validation.get("status") == "complete", "terminal manifest validation is incomplete", errors)
    require(validation.get("schema") == "valid", "terminal manifest schema is invalid", errors)
    require(validation.get("semantic") == "valid", "terminal manifest semantics are invalid", errors)
    require(validation.get("filesystem") == "tracked_inputs_verified", "terminal manifest filesystem validation failed", errors)
    artifact_types = Counter(artifact["artifact_type"] for artifact in terminal_manifest.get("artifacts", []))
    require(artifact_types["transcript_abundance"] == EXPECTED_SAMPLES, "terminal manifest lacks 23 transcript abundances", errors)
    require(artifact_types["differential_expression"] == EXPECTED_CONTRASTS, "terminal manifest lacks six contrasts", errors)

    required_html = [results / "execution_report.html", results / "timeline.html", results / "dag.html", multiqc]
    require(all(path.is_file() and path.stat().st_size > 0 for path in required_html), "required HTML report is missing", errors)

    trace_rows = read_tsv(results / "trace.tsv")
    statuses = Counter(row["status"] for row in trace_rows)
    require(statuses.get("FAILED", 0) == 0, "final trace contains failed tasks", errors)
    require(statuses.get("CACHED", 0) == 815, "unexpected cached-task count", errors)
    require(statuses.get("COMPLETED", 0) == 12, "unexpected completed-task count", errors)

    performance: dict[str, dict[str, float | int | Counter[str]]] = defaultdict(
        lambda: {"tasks": 0, "cached": 0, "realtime": 0.0, "max_realtime": 0.0, "max_rss": 0, "statuses": Counter()}
    )
    for row in trace_rows:
        process = row["name"].split(" (")[0]
        item = performance[process]
        runtime = seconds(row.get("realtime", ""))
        rss = memory_bytes(row.get("peak_rss", ""))
        item["tasks"] = int(item["tasks"]) + 1
        item["cached"] = int(item["cached"]) + (1 if row["status"] == "CACHED" else 0)
        item["realtime"] = float(item["realtime"]) + runtime
        item["max_realtime"] = max(float(item["max_realtime"]), runtime)
        item["max_rss"] = max(int(item["max_rss"]), rss)
        statuses_for_process = item["statuses"]
        assert isinstance(statuses_for_process, Counter)
        statuses_for_process[row["status"]] += 1

    with (output / "performance.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("process", "tasks", "cached", "completed", "total_realtime_seconds", "max_realtime_seconds", "max_peak_rss_bytes"))
        for process in sorted(performance):
            item = performance[process]
            process_statuses = item["statuses"]
            assert isinstance(process_statuses, Counter)
            writer.writerow((process, item["tasks"], item["cached"], process_statuses.get("COMPLETED", 0),
                             f"{float(item['realtime']):.3f}", f"{float(item['max_realtime']):.3f}", item["max_rss"]))

    qc_summary = {
        "fastqc_html": fastqc_count,
        "mapping_percent_min": min(mapping),
        "mapping_percent_mean": sum(mapping) / len(mapping),
        "mapping_percent_max": max(mapping),
        "multiqc": "PASS",
        "salmon_samples": len(mapping),
        "salmon_targets": sorted(targets),
        "salmon_versions": sorted(salmon_versions),
    }
    (output / "qc_summary.json").write_text(json.dumps(qc_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    final = {
        "schema_version": "1.0",
        "study": "PRJEB14695",
        "status": "PASS_WITH_LIMITATIONS" if not errors else "FAIL",
        "ready_for_review": not errors,
        "helixforge_commit": args.helixforge_commit,
        "runs": EXPECTED_RUNS,
        "samples": EXPECTED_SAMPLES,
        "fastq_files": EXPECTED_FASTQS,
        "matrices": matrices,
        "contrasts": contrast_rows,
        "terminal_artifact_types": dict(sorted(artifact_types.items())),
        "trace": dict(sorted(statuses.items())),
        "gates": {
            "fastq_integrity": "PASS" if integrity.get("status") == "PASS" else "FAIL",
            "technical_run_aggregation": "PASS" if accounting.get("status") == "PASS" else "FAIL",
            "qc": "PASS" if fastqc_count == EXPECTED_FASTQC_HTML and multiqc.is_file() else "FAIL",
            "quantification": "PASS" if len(mapping) == EXPECTED_SAMPLES else "FAIL",
            "expression_import": "PASS" if matrices["counts"] == {"genes": EXPECTED_GENES, "samples": EXPECTED_SAMPLES} else "FAIL",
            "differential_expression": "PASS" if len(contrast_rows) == EXPECTED_CONTRASTS else "FAIL",
            "terminal_manifest": "PASS" if validation.get("status") == "complete" else "FAIL",
            "html_reports": "PASS" if all(path.is_file() for path in required_html) else "FAIL",
            "resume_cache_reuse": "PASS" if statuses.get("CACHED") == 815 else "FAIL",
        },
        "limitations": [
            "The first full execution stopped at Import sample-table validation because technical_unit was incorrectly treated as invariant.",
            "HelixForge PR 88 corrected the run-level field contract; the controlled resume reused all 815 eligible upstream tasks.",
            "Candidate-gene reporting was not run because no candidate list was preregistered for this study.",
        ],
        "errors": errors,
    }
    (output / "final_validation.json").write_text(json.dumps(final, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": final["status"], "errors": errors}, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
