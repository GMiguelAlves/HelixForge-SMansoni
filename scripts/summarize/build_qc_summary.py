#!/usr/bin/env python3
"""Build the public, study-independent RNA-seq QC summary contract."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import zipfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Iterable


SAMPLE_COLUMNS = (
    "sample_id",
    "raw_reads_both_mates",
    "trimmed_reads_both_mates",
    "trim_retention_percent",
    "salmon_mapping_percent",
    "classification",
    "reason",
)


@dataclass(frozen=True)
class SourceFile:
    name: str
    text: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-id", required=True)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--results-root", type=Path)
    source.add_argument("--audit-zip", type=Path)
    parser.add_argument(
        "--sample-table",
        type=Path,
        help="Existing sample QC TSV to normalize; source results remain optional metadata evidence.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mapping-review-threshold", type=float, default=60.0)
    parser.add_argument("--retention-review-threshold", type=float, default=80.0)
    args = parser.parse_args()
    if not args.sample_table and not (args.results_root or args.audit_zip):
        parser.error("one of --sample-table, --results-root, or --audit-zip is required")
    return args


def read_source_files(results_root: Path | None, audit_zip: Path | None) -> list[SourceFile]:
    files: list[SourceFile] = []
    if results_root:
        for path in sorted(results_root.rglob("*")):
            if not path.is_file():
                continue
            posix = path.relative_to(results_root).as_posix()
            if posix.endswith("meta_info.json") or posix.endswith("_fastqc.html"):
                files.append(SourceFile(posix, path.read_text(encoding="utf-8")))
            elif "/multiqc/" in f"/{posix}" and posix.endswith(".html"):
                files.append(SourceFile(posix, ""))
        return files
    if audit_zip:
        with zipfile.ZipFile(audit_zip) as archive:
            for member in sorted(archive.namelist()):
                if member.endswith("meta_info.json") or member.endswith("_fastqc.html"):
                    files.append(SourceFile(member, archive.read(member).decode("utf-8")))
                elif "/multiqc/" in member and member.endswith(".html"):
                    files.append(SourceFile(member, ""))
    return files


def total_sequences(html: str, source: str) -> int:
    match = re.search(r"<td>Total Sequences</td><td>([0-9,]+)</td>", html)
    if not match:
        raise ValueError(f"FastQC Total Sequences is missing from {source}")
    return int(match.group(1).replace(",", ""))


def sample_from_meta_path(name: str) -> str:
    parent = PurePosixPath(name).parent.name
    sample = parent.split(".quantification", 1)[0]
    if "." in sample:
        sample = sample.split(".", 1)[1]
    if not sample:
        raise ValueError(f"cannot derive sample ID from {name}")
    return sample


def discover_metrics(files: Iterable[SourceFile]) -> tuple[dict[str, dict[str, float | int]], dict[str, object]]:
    file_list = list(files)
    salmon: dict[str, dict[str, float | int]] = {}
    versions: set[str] = set()
    targets: set[int] = set()
    for item in file_list:
        if not item.name.endswith("meta_info.json") or "/native_quantification/salmon_quant/" not in item.name:
            continue
        document = json.loads(item.text)
        sample_id = sample_from_meta_path(item.name)
        if sample_id in salmon:
            raise ValueError(f"duplicate Salmon metrics for {sample_id}")
        salmon[sample_id] = {
            "processed": int(document["num_processed"]),
            "mapped": int(document["num_mapped"]),
            "mapping": float(document["percent_mapped"]),
        }
        versions.add(str(document["salmon_version"]))
        targets.add(int(document["num_valid_targets"]))

    fastqc = {
        PurePosixPath(item.name).name: item
        for item in file_list
        if item.name.endswith("_fastqc.html") and "/native_qc/fastqc/" in item.name
    }
    sample_counts: dict[str, dict[str, int]] = {
        sample_id: {"raw": 0, "trimmed": 0} for sample_id in salmon
    }
    for sample_id in sorted(salmon):
        pattern = re.compile(
            rf"^{re.escape(sample_id)}_([A-Z]RR[0-9]+)_R([12])_trimmed_fastqc\.html$"
        )
        run_mates: list[tuple[str, str, SourceFile]] = []
        for basename, item in fastqc.items():
            match = pattern.fullmatch(basename)
            if match:
                run_mates.append((match.group(1), match.group(2), item))
        if not run_mates:
            raise ValueError(f"no run-level trimmed FastQC reports found for {sample_id}")
        seen: set[tuple[str, str]] = set()
        for run, mate, trimmed_item in run_mates:
            key = (run, mate)
            if key in seen:
                raise ValueError(f"duplicate trimmed FastQC report for {sample_id}/{run}/R{mate}")
            seen.add(key)
            raw_candidates = (
                f"{run}_{mate}_fastqc.html",
                f"{sample_id}_{run}_R{mate}_fastqc.html",
            )
            raw_item = next((fastqc[name] for name in raw_candidates if name in fastqc), None)
            if raw_item is None:
                raise ValueError(f"raw FastQC report missing for {sample_id}/{run}/R{mate}")
            sample_counts[sample_id]["raw"] += total_sequences(raw_item.text, raw_item.name)
            sample_counts[sample_id]["trimmed"] += total_sequences(trimmed_item.text, trimmed_item.name)

    metadata = {
        "fastqc_html": len(fastqc),
        "multiqc": "PASS" if any("/native_qc/multiqc/" in item.name and item.name.endswith(".html") for item in file_list) else "NOT_AVAILABLE",
        "salmon_targets": sorted(targets),
        "salmon_versions": sorted(versions),
        "processed_fragments_total": sum(int(row["processed"]) for row in salmon.values()),
        "mapped_fragments_total": sum(int(row["mapped"]) for row in salmon.values()),
    }
    rows: dict[str, dict[str, float | int]] = {}
    for sample_id in sorted(salmon):
        rows[sample_id] = {
            **sample_counts[sample_id],
            **salmon[sample_id],
        }
    return rows, metadata


def first_value(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = row.get(name, "").strip()
        if value:
            return value
    return ""


def load_sample_table(path: Path) -> dict[str, dict[str, float | int]]:
    with path.open(encoding="utf-8", newline="") as handle:
        table = list(csv.DictReader(handle, delimiter="\t"))
    rows: dict[str, dict[str, float | int]] = {}
    for row in table:
        sample_id = first_value(row, "sample_id")
        if not sample_id or sample_id in rows:
            raise ValueError(f"missing or duplicate sample_id in {path}: {sample_id!r}")
        raw_reads = first_value(row, "raw_reads_both_mates")
        trimmed_reads = first_value(row, "trimmed_reads_both_mates")
        if not raw_reads:
            raw_pairs = first_value(row, "raw_read_pairs")
            raw_reads = str(int(float(raw_pairs)) * 2) if raw_pairs else ""
        if not trimmed_reads:
            trimmed_pairs = first_value(row, "trimmed_read_pairs")
            trimmed_reads = str(int(float(trimmed_pairs)) * 2) if trimmed_pairs else ""
        mapping = first_value(row, "salmon_mapping_percent")
        if not raw_reads or not trimmed_reads or not mapping:
            raise ValueError(f"incomplete QC metrics for {sample_id} in {path}")
        raw = int(float(raw_reads))
        trimmed = int(float(trimmed_reads))
        rows[sample_id] = {
            "raw": raw,
            "trimmed": trimmed,
            "mapping": float(mapping),
            "processed": int(float(first_value(row, "salmon_processed_fragments", "salmon_fragments_processed") or 0)),
            "mapped": int(float(first_value(row, "salmon_mapped_fragments", "salmon_fragments_mapped") or 0)),
        }
    return rows


def classify(
    retention: float,
    mapping: float,
    retention_threshold: float,
    mapping_threshold: float,
) -> tuple[str, str]:
    reasons: list[str] = []
    if retention < retention_threshold:
        reasons.append("trim_retention_below_review_threshold")
    if mapping < mapping_threshold:
        reasons.append("salmon_mapping_below_review_threshold")
    return ("REVIEW", ";".join(reasons)) if reasons else ("PASS", "none")


def main() -> int:
    args = parse_args()
    source_files = read_source_files(args.results_root, args.audit_zip)
    discovered: dict[str, dict[str, float | int]] = {}
    metadata: dict[str, object] = {
        "fastqc_html": 0,
        "multiqc": "NOT_AVAILABLE",
        "salmon_targets": [],
        "salmon_versions": [],
        "processed_fragments_total": 0,
        "mapped_fragments_total": 0,
    }
    if source_files:
        discovered, metadata = discover_metrics(source_files)
    rows = load_sample_table(args.sample_table) if args.sample_table else discovered
    if not rows:
        raise ValueError("no sample QC metrics were discovered")
    if discovered:
        if set(rows) != set(discovered):
            raise ValueError("sample identities differ between sample table and preserved native results")
        for sample_id in rows:
            if not math.isclose(float(rows[sample_id]["mapping"]), float(discovered[sample_id]["mapping"]), abs_tol=1e-3):
                raise ValueError(f"Salmon mapping mismatch for {sample_id}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    sample_path = args.output_dir / "qc_sample_summary.tsv"
    classification_counts: dict[str, int] = {}
    retentions: list[float] = []
    mappings: list[float] = []
    with sample_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SAMPLE_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sample_id in sorted(rows):
            row = rows[sample_id]
            raw = int(row["raw"])
            trimmed = int(row["trimmed"])
            mapping = float(row["mapping"])
            if raw <= 0 or trimmed < 0 or not 0 <= mapping <= 100:
                raise ValueError(f"invalid QC metric for {sample_id}")
            retention = trimmed / raw * 100
            status, reason = classify(
                retention,
                mapping,
                args.retention_review_threshold,
                args.mapping_review_threshold,
            )
            retentions.append(retention)
            mappings.append(mapping)
            classification_counts[status] = classification_counts.get(status, 0) + 1
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "raw_reads_both_mates": raw,
                    "trimmed_reads_both_mates": trimmed,
                    "trim_retention_percent": f"{retention:.6f}",
                    "salmon_mapping_percent": f"{mapping:.6f}",
                    "classification": status,
                    "reason": reason,
                }
            )

    summary = {
        "schema_version": "1.0",
        "study_id": args.study_id,
        "samples": len(rows),
        "fastqc_html": int(metadata["fastqc_html"]),
        "multiqc": metadata["multiqc"],
        "raw_reads_total_both_mates": sum(int(row["raw"]) for row in rows.values()),
        "trimmed_reads_total_both_mates": sum(int(row["trimmed"]) for row in rows.values()),
        "trim_retention_percent_min": min(retentions),
        "trim_retention_percent_mean": statistics.mean(retentions),
        "trim_retention_percent_median": statistics.median(retentions),
        "trim_retention_percent_max": max(retentions),
        "processed_fragments_total": int(metadata["processed_fragments_total"]),
        "mapped_fragments_total": int(metadata["mapped_fragments_total"]),
        "salmon_samples": len(rows),
        "salmon_mapping_percent_min": min(mappings),
        "salmon_mapping_percent_mean": statistics.mean(mappings),
        "salmon_mapping_percent_median": statistics.median(mappings),
        "salmon_mapping_percent_max": max(mappings),
        "salmon_targets": metadata["salmon_targets"],
        "salmon_versions": metadata["salmon_versions"],
        "sample_classification": dict(sorted(classification_counts.items())),
        "thresholds": {
            "trim_retention_review_below_percent": args.retention_review_threshold,
            "salmon_mapping_review_below_percent": args.mapping_review_threshold,
        },
    }
    (args.output_dir / "qc_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"study_id": args.study_id, "samples": len(rows), "status": "complete"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
