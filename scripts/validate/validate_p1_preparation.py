#!/usr/bin/env python3
"""Validate the prepared P1 studies without requiring the raw FASTQs."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path


EXPECTED = {
    "PRJNA312093": (20, 20, "SINGLE", "full", "ENABLED"),
    "PRJNA494090": (24, 24, "PAIRED", "full", "ENABLED"),
    "PRJNA704646": (24, 8, "PAIRED", "import", "BLOCKED"),
    "PRJNA1150591": (15, 15, "PAIRED", "full", "ENABLED"),
    "PRJNA1177852": (8, 8, "PAIRED", "import", "BLOCKED"),
    "PRJNA551487": (10, 10, "PAIRED", "full", "ENABLED"),
    "PRJNA428470": (6, 6, "PAIRED", "full", "ENABLED"),
}
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def read(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate(root: Path) -> list[str]:
    errors: list[str] = []
    global_runs: set[str] = set()
    pending = {row["study"]: row for row in read(root / "metadata" / "p1_pending_manifest.tsv", "\t")}
    if set(pending) != set(EXPECTED):
        errors.append("p1_pending_manifest.tsv does not enumerate exactly the pending P1 studies")

    for study, (expected_runs, expected_samples, layout, mode, de_status) in EXPECTED.items():
        metadata_dir = root / "metadata" / study
        config_dir = root / "config" / study
        provenance_path = root / "provenance" / study / "input_manifest.json"
        snapshot = root / "provenance" / study / "source_snapshots" / "ena_read_run.tsv"
        required = [metadata_dir / "runs.tsv", metadata_dir / "samples.tsv", metadata_dir / "metadata.csv", provenance_path, snapshot, config_dir / "pipeline_config.sh", config_dir / "user_settings.sh", config_dir / "README.md"]
        missing = [str(path.relative_to(root)) for path in required if not path.is_file()]
        if missing:
            errors.append(f"{study}: missing files: {', '.join(missing)}")
            continue
        runs = read(metadata_dir / "runs.tsv", "\t")
        samples = read(metadata_dir / "samples.tsv", "\t")
        metadata = read(metadata_dir / "metadata.csv", ",")
        if (len(runs), len(samples), len(metadata)) != (expected_runs, expected_samples, expected_runs):
            errors.append(f"{study}: unexpected run/sample/metadata counts")
        run_ids = {row["run_accession"] for row in runs}
        if len(run_ids) != expected_runs:
            errors.append(f"{study}: duplicated run accessions")
        if global_runs & run_ids:
            errors.append(f"{study}: run accession reused across studies")
        global_runs.update(run_ids)
        if run_ids != {row["run_accession"] for row in metadata}:
            errors.append(f"{study}: runs.tsv and metadata.csv disagree")
        if {row["sample_id"] for row in runs} != {row["sample_id"] for row in samples}:
            errors.append(f"{study}: run/sample identity mismatch")
        counts = Counter(row["sample_id"] for row in runs)
        for sample in samples:
            if not SAFE_ID.fullmatch(sample["sample_id"]):
                errors.append(f"{study}: unsafe sample ID {sample['sample_id']}")
            if counts[sample["sample_id"]] != int(sample["technical_runs"]):
                errors.append(f"{study}: incorrect technical run count for {sample['sample_id']}")
        if {row["library_layout"] for row in runs} != {layout}:
            errors.append(f"{study}: expected homogeneous {layout} layout")
        for row in runs:
            fields = (row["fastq_1_url"], row["fastq_1_md5"], row["fastq_1_bytes"])
            if not all(fields) or not re.fullmatch(r"[0-9a-f]{32}", row["fastq_1_md5"]):
                errors.append(f"{study}/{row['run_accession']}: invalid R1 download contract")
            r2 = (row["fastq_2_url"], row["fastq_2_md5"], row["fastq_2_bytes"])
            if layout == "PAIRED" and (not all(r2) or not re.fullmatch(r"[0-9a-f]{32}", row["fastq_2_md5"])):
                errors.append(f"{study}/{row['run_accession']}: invalid R2 download contract")
            if layout == "SINGLE" and any(r2):
                errors.append(f"{study}/{row['run_accession']}: single-end run declares R2")
        for row in metadata:
            if row["library_layout"] != layout.lower():
                errors.append(f"{study}: metadata layout disagrees with runs.tsv")
            fragments = (row["fragment_length_mean"], row["fragment_length_sd"])
            if layout == "SINGLE" and not all(value.isdigit() and int(value) > 0 for value in fragments):
                errors.append(f"{study}: invalid single-end fragment prior")
            if layout == "PAIRED" and any(fragments):
                errors.append(f"{study}: paired record declares fragment prior")

        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        if provenance["source_snapshot"]["sha256"] != sha256(snapshot):
            errors.append(f"{study}: frozen ENA snapshot hash mismatch")
        if provenance["run_mode"] != mode or provenance["differential_expression"] != de_status:
            errors.append(f"{study}: execution decision differs from expected contract")
        settings = (config_dir / "user_settings.sh").read_text(encoding="utf-8")
        if f'export HF_RNASEQ_RUN_MODE="{mode}"' not in settings:
            errors.append(f"{study}: user_settings does not freeze run mode {mode}")

        if de_status == "BLOCKED":
            blocked_path = config_dir / "de_spec.blocked.json"
            if not blocked_path.is_file() or json.loads(blocked_path.read_text(encoding="utf-8")).get("status") != "BLOCKED":
                errors.append(f"{study}: explicit DE block is absent")
        else:
            spec_path = config_dir / "de_spec.json"
            if not spec_path.is_file():
                errors.append(f"{study}: DE specification is absent")
            else:
                spec = json.loads(spec_path.read_text(encoding="utf-8"))
                levels = Counter(row["condition"] for row in samples)
                for contrast in spec.get("contrasts", []):
                    for level in (contrast["numerator"], contrast["denominator"]):
                        if levels[level] < int(spec["parameters"]["min_replicates"]):
                            errors.append(f"{study}: contrast level {level} is under-replicated")

    if (root / "results" / "atlas").is_dir():
        # Existing historical atlas files are allowed, but this delivery must not
        # create a study-specific pending P1 integration tree.
        forbidden = [path for study in EXPECTED for path in (root / "results" / "atlas").rglob(f"*{study}*")]
        if forbidden:
            errors.append("pending P1 studies were written into results/atlas")
    return errors


def main() -> int:
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[2]
    errors = validate(root)
    print(json.dumps({"status": "PASS" if not errors else "FAIL", "errors": errors}, indent=2, ensure_ascii=False))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
