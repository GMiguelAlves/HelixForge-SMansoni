#!/usr/bin/env python3
"""Validate a registered study's downloaded paired-end FASTQs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


def digests(path: Path) -> tuple[str, str]:
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
    return md5.hexdigest(), sha256.hexdigest()


def validate_gzip(path: Path) -> None:
    with gzip.open(path, "rb") as handle:
        for _ in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            pass


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project")
    parser.add_argument("--package-root", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument("--output-sha256", required=True, type=Path)
    args = parser.parse_args()

    manifest = args.package_root / "metadata" / args.project / "runs.tsv"
    raw_dir = args.data_root / args.project / "fastq_ftp"
    with manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    errors: list[str] = []
    results: list[dict[str, object]] = []
    expected_names: set[str] = set()
    accessions = [row["run_accession"] for row in rows]
    if len(accessions) != len(set(accessions)):
        errors.append("duplicate run accessions in manifest")

    for row in rows:
        run = row["run_accession"]
        if row["library_layout"] != "PAIRED":
            errors.append(f"{run}: expected PAIRED layout")
        for mate in ("1", "2"):
            filename = f"{run}_{mate}.fastq.gz"
            expected_names.add(filename)
            path = raw_dir / filename
            record: dict[str, object] = {
                "run_accession": run,
                "sample_id": row["sample_id"],
                "mate": mate,
                "filename": filename,
                "expected_bytes": int(row[f"fastq_{mate}_bytes"]),
                "expected_md5": row[f"fastq_{mate}_md5"],
                "status": "FAIL",
            }
            if not path.is_file() or path.stat().st_size == 0:
                errors.append(f"{filename}: missing or empty")
                results.append(record)
                continue

            observed_bytes = path.stat().st_size
            observed_md5, observed_sha256 = digests(path)
            gzip_status = "PASS"
            try:
                validate_gzip(path)
            except (OSError, EOFError) as error:
                gzip_status = "FAIL"
                errors.append(f"{filename}: gzip integrity failed: {error}")

            record.update(
                {
                    "observed_bytes": observed_bytes,
                    "observed_md5": observed_md5,
                    "observed_sha256": observed_sha256,
                    "gzip": gzip_status,
                }
            )
            if observed_bytes != record["expected_bytes"]:
                errors.append(f"{filename}: size mismatch")
            if observed_md5 != record["expected_md5"]:
                errors.append(f"{filename}: MD5 mismatch")
            if (
                observed_bytes == record["expected_bytes"]
                and observed_md5 == record["expected_md5"]
                and gzip_status == "PASS"
            ):
                record["status"] = "PASS"
            results.append(record)

    observed_names = {path.name for path in raw_dir.glob("*.fastq.gz")}
    unexpected = sorted(observed_names - expected_names)
    missing = sorted(expected_names - observed_names)
    if unexpected:
        errors.append(f"unexpected FASTQs: {unexpected}")
    if missing:
        errors.append(f"missing FASTQs: {missing}")

    r1_bytes = sum(int(item.get("observed_bytes", 0)) for item in results if item["mate"] == "1")
    r2_bytes = sum(int(item.get("observed_bytes", 0)) for item in results if item["mate"] == "2")
    sizes = [int(item.get("observed_bytes", 0)) for item in results]
    report = {
        "schema_version": "1.0",
        "project": args.project,
        "status": "PASS" if not errors else "FAIL",
        "runs": len(rows),
        "samples": len({row["sample_id"] for row in rows}),
        "layout": "PAIRED",
        "fastq_files_expected": len(expected_names),
        "fastq_files_observed": len(observed_names),
        "pairing": "PASS" if not missing and len(expected_names) == 2 * len(rows) else "FAIL",
        "r1_total_bytes": r1_bytes,
        "r2_total_bytes": r2_bytes,
        "compressed_total_bytes": r1_bytes + r2_bytes,
        "largest_file_bytes": max(sizes, default=0),
        "smallest_file_bytes": min(sizes, default=0),
        "unexpected_files": unexpected,
        "missing_files": missing,
        "errors": errors,
        "files": results,
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    with args.output_tsv.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "run_accession", "sample_id", "mate", "filename", "expected_bytes",
            "observed_bytes", "expected_md5", "observed_md5", "observed_sha256",
            "gzip", "status",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    with args.output_sha256.open("w", encoding="utf-8") as handle:
        for item in sorted(results, key=lambda value: str(value["filename"])):
            if item.get("observed_sha256"):
                handle.write(f"{item['observed_sha256']}  {item['filename']}\n")

    print(json.dumps({key: value for key, value in report.items() if key != "files"}, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
