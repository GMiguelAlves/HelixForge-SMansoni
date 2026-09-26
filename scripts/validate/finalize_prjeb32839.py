#!/usr/bin/env python3
"""Validate and summarize the completed PRJEB32839 HelixForge run."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import zipfile
from collections import Counter, defaultdict
from pathlib import Path


EXPECTED_SAMPLES = 75
EXPECTED_RUNS = 150
EXPECTED_FASTQS = 300
EXPECTED_GENES = 9914
EXPECTED_TESTED_GENES = 9792
EXPECTED_CONTRASTS = 12
EXPECTED_FASTQC_HTML = 750
EXPECTED_TRACE_TASKS = 1153
EXPECTED_REPORT_QUERIES = 65
EXPECTED_REPORT_MATCHED = 48
EXPECTED_REPORT_UNMATCHED = 17


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


def fastqc_total_sequences(path: Path) -> int:
    with zipfile.ZipFile(path) as archive:
        member = next(name for name in archive.namelist() if name.endswith("/fastqc_data.txt"))
        with archive.open(member) as handle:
            for raw_line in handle:
                line = raw_line.decode("utf-8", errors="replace").rstrip("\n")
                if line.startswith("Total Sequences\t"):
                    return int(line.split("\t", 1)[1])
    raise ValueError(f"Total Sequences not found in {path.name}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--fastq-integrity", required=True, type=Path)
    parser.add_argument("--preflight", required=True, type=Path)
    parser.add_argument("--report-results", required=True, type=Path)
    parser.add_argument("--raw-fastqc-dir", required=True, type=Path)
    parser.add_argument("--trimmed-fastqc-dir", required=True, type=Path)
    parser.add_argument("--storage-summary", type=Path)
    parser.add_argument("--helixforge-commit", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    results = args.results
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    errors: list[str] = []

    integrity = json.loads(args.fastq_integrity.read_text(encoding="utf-8"))
    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    require(integrity.get("status") == "PASS", "FASTQ integrity did not pass", errors)
    require(integrity.get("runs") == EXPECTED_RUNS, "unexpected run count", errors)
    require(integrity.get("fastq_files_observed") == EXPECTED_FASTQS, "unexpected FASTQ count", errors)
    accounting = preflight.get("sample_accounting", {})
    require(preflight.get("status") == "PASS", "scientific preflight did not pass", errors)
    require(accounting.get("biological_samples") == EXPECTED_SAMPLES, "unexpected biological sample count", errors)
    require(accounting.get("deseq2_observations") == EXPECTED_SAMPLES, "DESeq2 observations are not biological samples", errors)
    require(accounting.get("technical_runs_are_independent_observations") is False,
            "technical runs were treated as independent observations", errors)

    import_root = results / "pipeline_info/native_import"
    matrix_root = import_root / "tximport"
    matrices: dict[str, dict[str, int]] = {}
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
    require(len(sample_rows) == EXPECTED_SAMPLES, "import sample table has an unexpected row count", errors)
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

    directions: dict[str, Counter[str]] = defaultdict(Counter)
    for row in read_tsv(de_root / "aggregate/differential_expression_results.tsv"):
        try:
            padj = float(row["padj"])
            log2fc = float(row["log2FoldChange"])
        except (TypeError, ValueError):
            continue
        if math.isfinite(padj) and padj < 0.05 and abs(log2fc) >= 1.0:
            directions[row["contrast"]]["up" if log2fc > 0 else "down"] += 1

    with (output / "differential_expression_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("contrast", "samples", "tested_genes", "significant", "up", "down", "criteria"))
        for row in contrast_rows:
            counts = directions[row["contrast"]]
            writer.writerow((row["contrast"], row["n_samples"], row["n_genes"], row["n_significant"],
                             counts["up"], counts["down"], "padj<0.05;abs(log2FoldChange)>=1"))
            require(counts["up"] + counts["down"] == int(row["n_significant"]),
                    f"significant-gene direction counts differ for {row['contrast']}", errors)

    quant_root = results / "pipeline_info/native_quantification/salmon_quant"
    meta_files = sorted(quant_root.glob("*.quantification.quantification_statistics/meta_info.json"))
    require(len(meta_files) == EXPECTED_SAMPLES, "missing Salmon meta_info records", errors)
    mapping: list[float] = []
    mapping_by_sample: dict[str, float] = {}
    processed_fragments: list[int] = []
    low_mapping: list[dict[str, float | int | str]] = []
    targets: set[int] = set()
    salmon_versions: set[str] = set()
    for path in meta_files:
        document = json.loads(path.read_text(encoding="utf-8"))
        percent = float(document["percent_mapped"])
        processed = int(document["num_processed"])
        sample = path.parent.name.split(".quantification", 1)[0]
        sample_id = sample.split(".", 1)[1] if "." in sample else sample
        mapping.append(percent)
        mapping_by_sample[sample_id] = percent
        processed_fragments.append(processed)
        if percent < 50:
            low_mapping.append({"sample": sample, "percent_mapped": percent, "processed_fragments": processed})
        targets.add(int(document["num_valid_targets"]))
        salmon_versions.add(str(document["salmon_version"]))
    require(mapping and all(math.isfinite(value) and 0 <= value <= 100 for value in mapping),
            "invalid Salmon mapping percentages", errors)
    require(targets == {10913}, f"unexpected Salmon target counts: {sorted(targets)}", errors)
    require(salmon_versions == {"1.10.3"}, f"unexpected Salmon versions: {sorted(salmon_versions)}", errors)

    fastqc_count = len(list((results / "pipeline_info/native_qc/fastqc").glob("*.html")))
    multiqc_candidates = list((results / "pipeline_info/native_qc/multiqc").glob("*.html"))
    require(fastqc_count == EXPECTED_FASTQC_HTML, f"unexpected FastQC HTML count: {fastqc_count}", errors)
    require(len(multiqc_candidates) == 1 and multiqc_candidates[0].stat().st_size > 0, "MultiQC HTML is missing", errors)

    raw_fastqc: dict[tuple[str, str], int] = {}
    for path in sorted(args.raw_fastqc_dir.glob("*.zip")):
        match = re.fullmatch(r"(ERR[0-9]+)_([12])_fastqc\.zip", path.name)
        if match:
            raw_fastqc[(match.group(1), match.group(2))] = fastqc_total_sequences(path)
    trimmed_fastqc: dict[tuple[str, str], int] = {}
    run_to_sample: dict[str, str] = {}
    for path in sorted(args.trimmed_fastqc_dir.glob("*.zip")):
        match = re.fullmatch(r"(.+)_(ERR[0-9]+)_R([12])_trimmed_fastqc\.zip", path.name)
        if match:
            sample_id, run, mate = match.groups()
            trimmed_fastqc[(run, mate)] = fastqc_total_sequences(path)
            previous = run_to_sample.setdefault(run, sample_id)
            require(previous == sample_id, f"run {run} maps to multiple samples in FastQC outputs", errors)
    require(len(raw_fastqc) == EXPECTED_FASTQS, f"unexpected raw FastQC data count: {len(raw_fastqc)}", errors)
    require(len(trimmed_fastqc) == EXPECTED_FASTQS,
            f"unexpected trimmed FastQC data count: {len(trimmed_fastqc)}", errors)
    require(set(raw_fastqc) == set(trimmed_fastqc), "raw/trimmed FastQC keys differ", errors)
    retention = [trimmed_fastqc[key] / raw_fastqc[key] * 100 for key in raw_fastqc if raw_fastqc[key] > 0]
    sample_reads: dict[str, Counter[str]] = defaultdict(Counter)
    for (run, mate), raw_total in raw_fastqc.items():
        sample_id = run_to_sample[run]
        sample_reads[sample_id]["raw"] += raw_total
        sample_reads[sample_id]["trimmed"] += trimmed_fastqc[(run, mate)]
    require(set(sample_reads) == set(mapping_by_sample), "FastQC and Salmon sample identities differ", errors)
    qc_classification = Counter()
    with (output / "qc_sample_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("sample_id", "raw_reads_both_mates", "trimmed_reads_both_mates",
                         "trim_retention_percent", "salmon_mapping_percent", "classification", "reason"))
        for sample_id in sorted(sample_reads):
            values = sample_reads[sample_id]
            sample_retention = values["trimmed"] / values["raw"] * 100
            sample_mapping = mapping_by_sample[sample_id]
            reasons = []
            if sample_retention < 80:
                reasons.append("trim_retention_below_80_percent")
            if sample_mapping < 50:
                reasons.append("salmon_mapping_below_50_percent")
            classification = "REVIEW" if reasons else "PASS"
            qc_classification[classification] += 1
            writer.writerow((sample_id, values["raw"], values["trimmed"], f"{sample_retention:.6f}",
                             f"{sample_mapping:.6f}", classification, ";".join(reasons)))

    validation = json.loads((results / "pipeline_info/integration_api/run_manifest.validation.json").read_text(encoding="utf-8"))
    terminal_manifest = json.loads((results / "pipeline_info/integration_api/rnaseq_run_manifest.json").read_text(encoding="utf-8"))
    require(validation.get("status") == "complete", "terminal manifest validation is incomplete", errors)
    require(validation.get("schema") == "valid", "terminal manifest schema is invalid", errors)
    require(validation.get("semantic") == "valid", "terminal manifest semantics are invalid", errors)
    require(validation.get("filesystem") == "tracked_inputs_verified", "terminal manifest filesystem validation failed", errors)
    artifact_types = Counter(artifact["artifact_type"] for artifact in terminal_manifest.get("artifacts", []))
    require(artifact_types["transcript_abundance"] == EXPECTED_SAMPLES,
            "terminal manifest lacks 75 transcript abundances", errors)
    require(artifact_types["differential_expression"] == EXPECTED_CONTRASTS,
            "terminal manifest lacks 12 differential-expression artifacts", errors)

    report = args.report_results
    report_context = json.loads((report / "context.json").read_text(encoding="utf-8"))
    report_execution = json.loads((report / "execution.json").read_text(encoding="utf-8"))
    report_manifest = json.loads((report / "manifest.json").read_text(encoding="utf-8"))
    catalog = read_tsv(report / "tables/gene_catalog.tsv")
    report_matched = sum(row["match_type"] != "unmatched" for row in catalog)
    report_unmatched = sum(row["match_type"] == "unmatched" for row in catalog)
    report_png = len(list(report.rglob("*.png")))
    require(report_context.get("status") == "complete", "gene-report context is incomplete", errors)
    require(report_execution.get("status") == "complete", "gene-report execution is incomplete", errors)
    require(report_context.get("sample_count") == EXPECTED_SAMPLES, "gene report has wrong sample count", errors)
    require(report_context.get("query_count") == EXPECTED_REPORT_QUERIES, "gene report has wrong query count", errors)
    require(report_matched == EXPECTED_REPORT_MATCHED, "gene report has wrong matched-gene count", errors)
    require(report_unmatched == EXPECTED_REPORT_UNMATCHED, "gene report has wrong unmatched-gene count", errors)
    require((report / "gene_set_report.html").stat().st_size > 0, "gene-report HTML is missing", errors)
    require(report_png == 348, f"unexpected gene-report figure count: {report_png}", errors)
    require(report_manifest.get("artifacts", {}).get("html", {}).get("available") is True,
            "gene-report manifest does not declare HTML", errors)

    required_html = [results / "execution_report.html", results / "timeline.html", results / "dag.html"] + multiqc_candidates
    require(all(path.is_file() and path.stat().st_size > 0 for path in required_html), "required HTML report is missing", errors)

    trace_rows = read_tsv(results / "trace.tsv")
    statuses = Counter(row["status"] for row in trace_rows)
    require(len(trace_rows) == EXPECTED_TRACE_TASKS, f"unexpected trace task count: {len(trace_rows)}", errors)
    require(statuses.get("FAILED", 0) == 0, "final trace contains failed tasks", errors)
    require(statuses.get("COMPLETED", 0) == EXPECTED_TRACE_TASKS, "not all scientific tasks completed", errors)

    performance: dict[str, dict[str, object]] = defaultdict(
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
        process_statuses = item["statuses"]
        assert isinstance(process_statuses, Counter)
        process_statuses[row["status"]] += 1

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
        "mapping_percent_mean": statistics.mean(mapping),
        "mapping_percent_median": statistics.median(mapping),
        "mapping_percent_max": max(mapping),
        "multiqc": "PASS",
        "raw_reads_total_both_mates": sum(raw_fastqc.values()),
        "trimmed_reads_total_both_mates": sum(trimmed_fastqc.values()),
        "trim_retention_percent_min": min(retention),
        "trim_retention_percent_mean": statistics.mean(retention),
        "trim_retention_percent_median": statistics.median(retention),
        "trim_retention_percent_max": max(retention),
        "processed_fragments_total": sum(processed_fragments),
        "salmon_samples": len(mapping),
        "salmon_targets": sorted(targets),
        "salmon_versions": sorted(salmon_versions),
        "samples_below_50_percent_mapping": low_mapping,
        "sample_classification": dict(sorted(qc_classification.items())),
    }
    (output / "qc_summary.json").write_text(json.dumps(qc_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    storage: dict[str, int] = {}
    if args.storage_summary and args.storage_summary.is_file():
        storage = {row["item"]: int(row["bytes"]) for row in read_tsv(args.storage_summary)}
        input_bytes = int(integrity["compressed_total_bytes"])
        storage["observed_pre_cleanup_multiplier"] = round(storage.get("project_total", 0) / input_bytes, 6)

    limitations = [
        "Fifteen 32-day sporocyst samples mapped below 50%; the stage-specific pattern is reported for biological review rather than silently excluded.",
        "Seventeen preregistered candidate identifiers were absent from the WBPS19 expression/annotation inputs and remain explicit unmatched entries.",
        "The candidate-gene report was produced with the official report_reentry mode after the full workflow because reporting was disabled in the original full invocation; no scientific process was recomputed.",
        "Storage was measured at checkpoints and immediately before cleanup, not sampled continuously; the observed pre-cleanup footprint is not claimed as an exact peak.",
    ]
    final = {
        "schema_version": "1.0",
        "study": "PRJEB32839",
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
        "gene_report": {
            "queries": len(catalog),
            "matched": report_matched,
            "unmatched": report_unmatched,
            "figures": report_png,
            "provider": report_context.get("provider"),
            "container": report_execution.get("container"),
        },
        "storage": storage,
        "gates": {
            "fastq_integrity": "PASS" if integrity.get("status") == "PASS" else "FAIL",
            "technical_run_aggregation": "PASS" if preflight.get("gates", {}).get("technical_run_aggregation") == "PASS" else "FAIL",
            "qc": "PASS" if fastqc_count == EXPECTED_FASTQC_HTML and len(multiqc_candidates) == 1 else "FAIL",
            "quantification": "PASS" if len(mapping) == EXPECTED_SAMPLES else "FAIL",
            "expression_import": "PASS" if matrices["counts"] == {"genes": EXPECTED_GENES, "samples": EXPECTED_SAMPLES} else "FAIL",
            "differential_expression": "PASS" if len(contrast_rows) == EXPECTED_CONTRASTS else "FAIL",
            "gene_report": "PASS" if report_execution.get("status") == "complete" and report_png == 348 else "FAIL",
            "terminal_manifest": "PASS" if validation.get("status") == "complete" else "FAIL",
            "html_reports": "PASS" if all(path.is_file() for path in required_html) else "FAIL",
            "resume_cache_reuse": "NOT_APPLICABLE_CLEAN_RUN",
        },
        "limitations": limitations,
        "errors": errors,
    }
    (output / "final_validation.json").write_text(json.dumps(final, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": final["status"], "errors": errors}, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
