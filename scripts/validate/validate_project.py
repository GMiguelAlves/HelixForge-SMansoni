#!/usr/bin/env python3
"""Validate the frozen HelixForge-SMansoni project contracts."""

from __future__ import annotations

import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path


EXPECTED = {
    "PRJNA602528": (10, 10),
    "PRJNA597909": (20, 20),
    "PRJEB14695": (138, 23),
    "PRJEB32839": (150, 75),
}
MASTER_FIELDS = [
    "study", "sample", "run", "layout", "stage", "sex", "tissue",
    "condition", "treatment", "time", "replicate", "batch", "source",
]


def read_rows(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def validate(root: Path) -> list[str]:
    errors: list[str] = []
    all_run_ids: set[str] = set()

    for project, (expected_runs, expected_samples) in EXPECTED.items():
        metadata_dir = root / "metadata" / project
        config_dir = root / "config" / project
        runs = read_rows(metadata_dir / "runs.tsv", "\t")
        samples = read_rows(metadata_dir / "samples.tsv", "\t")
        metadata = read_rows(metadata_dir / "metadata.csv", ",")

        if (len(runs), len(samples), len(metadata)) != (
            expected_runs,
            expected_samples,
            expected_runs,
        ):
            errors.append(f"{project}: unexpected run/sample/metadata counts")

        run_ids = {row["run_accession"] for row in runs}
        if len(run_ids) != expected_runs:
            errors.append(f"{project}: duplicate run accessions")
        if all_run_ids & run_ids:
            errors.append(f"{project}: run accession reused across studies")
        all_run_ids.update(run_ids)

        if run_ids != {row["run_accession"] for row in metadata}:
            errors.append(f"{project}: run metadata mismatch")
        if {row["sample_id"] for row in samples} != {row["sample_id"] for row in runs}:
            errors.append(f"{project}: sample identity mismatch")

        counts = Counter(row["sample_id"] for row in runs)
        for sample in samples:
            sample_id = sample["sample_id"]
            if counts[sample_id] != int(sample["technical_runs"]):
                errors.append(f"{project}: technical run count for {sample_id}")
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", sample_id):
                errors.append(f"{project}: unsafe sample identifier {sample_id}")

        if any(row["library_layout"] != "PAIRED" for row in runs):
            errors.append(f"{project}: non-paired run in the frozen input set")

        if project == "PRJNA602528":
            blocked = json.loads((config_dir / "de_spec.blocked.json").read_text(encoding="utf-8"))
            if blocked.get("status") != "BLOCKED":
                errors.append(f"{project}: explicit DE block missing")
        else:
            spec = json.loads((config_dir / "de_spec.json").read_text(encoding="utf-8"))
            levels = Counter(sample["condition"] for sample in samples)
            for contrast in spec["contrasts"]:
                for level in (contrast["numerator"], contrast["denominator"]):
                    if levels[level] < spec["parameters"]["min_replicates"]:
                        errors.append(f"{project}: under-replicated level {level}")

    master = read_rows(root / "metadata" / "master_sample_manifest.tsv", "\t")
    if not master or list(master[0]) != MASTER_FIELDS:
        errors.append("master manifest: unexpected columns")
    if len(master) != sum(value[0] for value in EXPECTED.values()):
        errors.append("master manifest: unexpected row count")
    if {row["run"] for row in master} != all_run_ids:
        errors.append("master manifest: run identities do not match study manifests")

    reference = json.loads(
        (root / "config" / "reference" / "reference_manifest.json").read_text(encoding="utf-8")
    )
    if reference.get("assembly") != "SM_V10":
        errors.append("reference manifest: assembly is not SM_V10")
    if reference.get("wormbase_parasite_release") != "WBPS19":
        errors.append("reference manifest: WormBase ParaSite release is not WBPS19")
    if {item.get("role") for item in reference.get("artifacts", [])} != {
        "genome", "transcriptome", "annotation_gff3", "annotation_gtf",
    }:
        errors.append("reference manifest: incomplete artifact roles")
    if reference.get("status") != "VERIFIED":
        errors.append("reference manifest: reference is not verified")
    if reference.get("checksums", {}).get("status") != "VERIFIED":
        errors.append("reference manifest: checksum set is not verified")
    for item in reference.get("artifacts", []):
        if item.get("checksum_status") != "VERIFIED":
            errors.append(f"reference manifest: unverified {item.get('role')} artifact")
        if not re.fullmatch(r"[0-9a-f]{64}", item.get("sha256", "")):
            errors.append(f"reference manifest: invalid {item.get('role')} SHA-256")

    reference_validation = json.loads(
        (root / "provenance" / "reference" / "reference_validation.json").read_text(
            encoding="utf-8"
        )
    )
    if reference_validation.get("status") != "PASS":
        errors.append("reference validation: status is not PASS")
    if reference_validation.get("contig_compatibility", {}).get("status") != "PASS":
        errors.append("reference validation: contig compatibility failed")

    server_state = json.loads(
        (root / "provenance" / "server_preparation_state.json").read_text(encoding="utf-8")
    )
    if server_state.get("prjna602528", {}).get("mode") != "IMPORT_ONLY":
        errors.append("server preparation: PRJNA602528 is not import-only")
    if server_state.get("prjna602528", {}).get("differential_expression") != "BLOCKED":
        errors.append("server preparation: differential expression is not blocked")

    return errors


def main() -> int:
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[2]
    errors = validate(root)
    print(json.dumps({"status": "PASS" if not errors else "FAIL", "errors": errors}, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
