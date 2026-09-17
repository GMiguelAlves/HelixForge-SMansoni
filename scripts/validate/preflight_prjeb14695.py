#!/usr/bin/env python3
"""Scientific preflight for PRJEB14695 run aggregation and DE design."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlparse


EXPECTED_RUNS = 138
EXPECTED_SAMPLES = 23
EXPECTED_RUNS_PER_SAMPLE = 6
EXPECTED_STUDY = "PRJEB14695"
EXPECTED_SECONDARY_STUDY = "ERP016356"
AGGREGATION_POLICY = "trim_each_run_then_concatenate_by_mate_before_salmon"
def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path: Path, delimiter: str = "\t") -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return list(reader.fieldnames or []), [dict(row) for row in reader]


def write_table(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def split_field(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def normalized_remote_path(value: str) -> str:
    parsed = urlparse(value if "://" in value else f"ftp://{value}")
    return f"{parsed.netloc.lower()}{parsed.path}"


def filename(value: str) -> str:
    return Path(urlparse(value if "://" in value else f"ftp://{value}").path).name


def matrix_rank(matrix: list[list[float]], tolerance: float = 1e-10) -> int:
    work = [row[:] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if work else 0
    rank = 0
    for column in range(columns):
        pivot = max(range(rank, rows), key=lambda row: abs(work[row][column]), default=rank)
        if rank >= rows or abs(work[pivot][column]) <= tolerance:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        work[rank] = [value / pivot_value for value in work[rank]]
        for row in range(rows):
            if row != rank and abs(work[row][column]) > tolerance:
                factor = work[row][column]
                work[row] = [a - factor * b for a, b in zip(work[row], work[rank])]
        rank += 1
    return rank


def add_error(errors: list[str], condition: bool, message: str) -> None:
    if condition:
        errors.append(message)


def validate(root: Path, ena_path: Path, output: Path) -> dict[str, object]:
    errors: list[str] = []
    metadata_dir = root / "metadata" / EXPECTED_STUDY
    _, runs = read_table(metadata_dir / "runs.tsv")
    _, samples = read_table(metadata_dir / "samples.tsv")
    _, metadata = read_table(metadata_dir / "metadata.csv", ",")
    _, ena = read_table(ena_path)
    spec = json.loads((root / "config" / EXPECTED_STUDY / "de_spec.json").read_text(encoding="utf-8"))

    add_error(errors, len(runs) != EXPECTED_RUNS, f"frozen run count is {len(runs)}, expected {EXPECTED_RUNS}")
    add_error(errors, len(metadata) != EXPECTED_RUNS, f"metadata row count is {len(metadata)}, expected {EXPECTED_RUNS}")
    add_error(errors, len(ena) != EXPECTED_RUNS, f"ENA run count is {len(ena)}, expected {EXPECTED_RUNS}")
    add_error(errors, len(samples) != EXPECTED_SAMPLES, f"sample count is {len(samples)}, expected {EXPECTED_SAMPLES}")

    run_counts = Counter(row["run_accession"] for row in runs)
    metadata_counts = Counter(row["run_accession"] for row in metadata)
    ena_counts = Counter(row["run_accession"] for row in ena)
    duplicated_runs = sorted(run for run, count in run_counts.items() if count != 1)
    duplicated_metadata = sorted(run for run, count in metadata_counts.items() if count != 1)
    duplicated_ena = sorted(run for run, count in ena_counts.items() if count != 1)
    add_error(errors, bool(duplicated_runs), f"frozen duplicate run accessions: {duplicated_runs[:10]}")
    add_error(errors, bool(duplicated_metadata), f"metadata duplicate run accessions: {duplicated_metadata[:10]}")
    add_error(errors, bool(duplicated_ena), f"ENA duplicate run accessions: {duplicated_ena[:10]}")

    frozen_by_run = {row["run_accession"]: row for row in runs}
    metadata_by_run = {row["run_accession"]: row for row in metadata}
    ena_by_run = {row["run_accession"]: row for row in ena}
    frozen_set, metadata_set, ena_set = set(frozen_by_run), set(metadata_by_run), set(ena_by_run)
    add_error(errors, frozen_set != metadata_set, "frozen runs and workflow metadata contain different run accessions")
    add_error(errors, frozen_set != ena_set, "frozen runs and live ENA contain different run accessions")
    ena_queries = {row.get("ena_query_url", "") for row in runs}
    add_error(errors, len(ena_queries) != 1 or not next(iter(ena_queries), ""), "frozen runs do not share one ENA query URL")

    sample_by_id = {row["sample_id"]: row for row in samples}
    add_error(errors, len(sample_by_id) != len(samples), "sample table contains duplicate sample_id values")
    run_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    conflicts: list[str] = []
    ena_differences: list[str] = []
    for accession in sorted(frozen_set):
        frozen = frozen_by_run[accession]
        workflow = metadata_by_run.get(accession)
        public = ena_by_run.get(accession)
        sample_id = frozen["sample_id"]
        run_groups[sample_id].append(frozen)
        if sample_id not in sample_by_id:
            conflicts.append(f"{accession}: unknown sample_id {sample_id}")
            continue
        sample = sample_by_id[sample_id]
        if workflow is None or public is None:
            continue
        for field in ("sample_id", "biosample", "library_name", "technical_unit"):
            if workflow.get(field) != frozen.get(field):
                conflicts.append(f"{accession}: {field} differs between runs.tsv and metadata.csv")
        for field in ("condition", "stage", "tissue", "sex", "infection_mode", "treatment", "time_hours", "batch", "replicate"):
            if workflow.get(field) != sample.get(field):
                conflicts.append(f"{accession}: {field} differs from biological sample {sample_id}")

        checks = {
            "study_accession": EXPECTED_STUDY,
            "secondary_study_accession": EXPECTED_SECONDARY_STUDY,
            "experiment_accession": frozen["experiment_accession"],
            "sample_accession": frozen["biosample"],
            "secondary_sample_accession": frozen["secondary_biosample"],
            "sample_alias": sample["source_sample_name"],
            "library_name": frozen["library_name"],
            "library_layout": "PAIRED",
            "library_strategy": "RNA-Seq",
            "library_source": "TRANSCRIPTOMIC",
        }
        for field, expected in checks.items():
            if public.get(field, "") != expected:
                ena_differences.append(f"{accession}: ENA {field}={public.get(field, '')!r}, frozen={expected!r}")
        for field, frozen_fields in (
            ("fastq_ftp", ("fastq_1_url", "fastq_2_url")),
            ("fastq_md5", ("fastq_1_md5", "fastq_2_md5")),
            ("fastq_bytes", ("fastq_1_bytes", "fastq_2_bytes")),
        ):
            public_values = split_field(public.get(field, ""))
            expected_values = [frozen[name] for name in frozen_fields]
            if field == "fastq_ftp":
                public_values = [normalized_remote_path(value) for value in public_values]
                expected_values = [normalized_remote_path(value) for value in expected_values]
            if public_values != expected_values:
                ena_differences.append(f"{accession}: ENA {field} differs from frozen inventory")
        for field in ("read_count", "base_count", "first_public"):
            frozen_field = "ena_first_public" if field == "first_public" else field
            if public.get(field, "") != frozen.get(frozen_field, ""):
                ena_differences.append(f"{accession}: ENA {field} differs from frozen inventory")

    add_error(errors, bool(conflicts), "run-to-sample conflicts: " + "; ".join(conflicts[:10]))
    add_error(errors, bool(ena_differences), "live ENA differences: " + "; ".join(ena_differences[:10]))
    add_error(errors, set(run_groups) != set(sample_by_id), "run inventory and sample inventory contain different sample IDs")

    aggregation_rows: list[dict[str, object]] = []
    mapping_rows: list[dict[str, object]] = []
    for sample_id in sorted(sample_by_id):
        sample = sample_by_id[sample_id]
        grouped = sorted(run_groups.get(sample_id, []), key=lambda row: row["run_accession"])
        technical_units = [row["technical_unit"] for row in grouped]
        expected_units = [f"technical_run_{index}" for index in range(1, EXPECTED_RUNS_PER_SAMPLE + 1)]
        if len(grouped) != EXPECTED_RUNS_PER_SAMPLE:
            errors.append(f"{sample_id}: has {len(grouped)} runs, expected {EXPECTED_RUNS_PER_SAMPLE}")
        if sorted(technical_units) != expected_units:
            errors.append(f"{sample_id}: invalid or duplicate technical-unit inventory")
        if int(sample["technical_runs"]) != len(grouped):
            errors.append(f"{sample_id}: samples.tsv technical_runs disagrees with mapped runs")
        for row in grouped:
            mapping_rows.append({
                "run_accession": row["run_accession"],
                "biological_sample_id": sample_id,
                "biosample": row["biosample"],
                "technical_unit": row["technical_unit"],
                "r1_file": filename(row["fastq_1_url"]),
                "r2_file": filename(row["fastq_2_url"]),
                "condition": sample["condition"],
                "sex": sample["sex"],
                "tissue_context": sample["tissue"],
                "pairing_state": sample["infection_mode"],
                "replicate": sample["replicate"],
            })
        aggregation_rows.append({
            "biological_sample_id": sample_id,
            "number_of_runs": len(grouped),
            "run_accessions": ";".join(row["run_accession"] for row in grouped),
            "r1_files": ";".join(filename(row["fastq_1_url"]) for row in grouped),
            "r2_files": ";".join(filename(row["fastq_2_url"]) for row in grouped),
            "biological_group": sample["condition"],
            "sex": sample["sex"],
            "tissue_context": sample["tissue"],
            "pairing_state": sample["infection_mode"],
            "replicate": sample["replicate"],
            "aggregation_policy": AGGREGATION_POLICY,
            "deseq2_observation": sample_id,
        })

    design = spec.get("design", {})
    variable = str(design.get("variable", ""))
    covariates = [str(value) for value in design.get("covariates", [])]
    formula = str(design.get("formula", ""))
    expected_formula = "~ " + " + ".join(covariates + [variable])
    add_error(errors, formula != expected_formula, f"design formula {formula!r} does not match declared fields")
    add_error(errors, variable not in (samples[0] if samples else {}), f"design variable {variable!r} is absent")
    levels = sorted({row.get(variable, "") for row in samples})
    add_error(errors, "" in levels, "design variable contains an empty level")
    reference_level = levels[0] if levels else ""
    design_columns = ["intercept"] + [f"{variable}[{level}]" for level in levels[1:]]
    design_rows: list[dict[str, object]] = []
    matrix: list[list[float]] = []
    for sample in sorted(samples, key=lambda row: row["sample_id"]):
        values = [1.0] + [1.0 if sample[variable] == level else 0.0 for level in levels[1:]]
        matrix.append(values)
        design_rows.append({"sample_id": sample["sample_id"], **dict(zip(design_columns, map(int, values)))})
    rank = matrix_rank(matrix)
    full_rank = bool(matrix) and rank == len(design_columns)
    add_error(errors, not full_rank, f"design matrix rank {rank} is below {len(design_columns)}")

    level_counts = Counter(row.get(variable, "") for row in samples)
    contrast_rows: list[dict[str, object]] = []
    contrast_ids: set[str] = set()
    min_replicates = int(spec.get("parameters", {}).get("min_replicates", 0))
    for contrast in spec.get("contrasts", []):
        contrast_id = str(contrast.get("id", ""))
        numerator = str(contrast.get("numerator", ""))
        denominator = str(contrast.get("denominator", ""))
        estimable = (
            contrast.get("factor") == variable
            and numerator in level_counts
            and denominator in level_counts
            and numerator != denominator
            and level_counts[numerator] >= min_replicates
            and level_counts[denominator] >= min_replicates
            and contrast_id not in contrast_ids
        )
        contrast_ids.add(contrast_id)
        if not estimable:
            errors.append(f"contrast {contrast_id!r} is not estimable under the frozen design")
        contrast_rows.append({
            "contrast": contrast_id,
            "factor": variable,
            "numerator": numerator,
            "denominator": denominator,
            "positive_log2fc_means": f"higher in {numerator}",
            "negative_log2fc_means": f"higher in {denominator}",
            "numerator_replicates": level_counts[numerator],
            "denominator_replicates": level_counts[denominator],
            "estimable": str(estimable).upper(),
        })

    accession_conflicts = []
    for path in [metadata_dir / "runs.tsv", metadata_dir / "samples.tsv", metadata_dir / "metadata.csv", root / "config" / EXPECTED_STUDY / "README.md", root / "config" / EXPECTED_STUDY / "de_spec.json"]:
        if "PRJEB1237" in path.read_text(encoding="utf-8"):
            accession_conflicts.append(str(path.relative_to(root)))
    add_error(errors, bool(accession_conflicts), f"forbidden PRJEB1237 reference in {accession_conflicts}")

    write_table(output / "run_to_sample_mapping.tsv", list(mapping_rows[0]), mapping_rows)
    write_table(output / "technical_run_aggregation.tsv", list(aggregation_rows[0]), aggregation_rows)
    write_table(output / "design_matrix.tsv", list(design_rows[0]), design_rows)
    write_table(output / "contrast_orientation.tsv", list(contrast_rows[0]), contrast_rows)

    gates = {
        "run_to_sample_mapping": "PASS" if not conflicts and not duplicated_runs and frozen_set == metadata_set == ena_set else "FAIL",
        "technical_run_aggregation": "PASS" if all(len(rows) == EXPECTED_RUNS_PER_SAMPLE for rows in run_groups.values()) and len(run_groups) == EXPECTED_SAMPLES else "FAIL",
        "biological_sample_accounting": "PASS" if len(aggregation_rows) == EXPECTED_SAMPLES and sum(int(row["number_of_runs"]) for row in aggregation_rows) == EXPECTED_RUNS else "FAIL",
        "ena_reconciliation": "PASS" if not ena_differences and frozen_set == ena_set else "FAIL",
        "accession_conflict": "PASS" if not accession_conflicts else "FAIL",
        "statistical_design_preflight": "PASS" if full_rank and all(row["estimable"] == "TRUE" for row in contrast_rows) else "FAIL",
    }
    if errors:
        for name in gates:
            if gates[name] == "PASS" and name in {"run_to_sample_mapping", "technical_run_aggregation", "biological_sample_accounting"}:
                gates[name] = "FAIL"
    report: dict[str, object] = {
        "schema_version": "1.0",
        "study": EXPECTED_STUDY,
        "status": "PASS" if not errors and all(value == "PASS" for value in gates.values()) else "FAIL",
        "gates": gates,
        "run_inventory": {
            "expected": EXPECTED_RUNS,
            "frozen": len(runs),
            "workflow_metadata": len(metadata),
            "ena_live": len(ena),
            "mapped": len(mapping_rows),
            "unmapped": len(frozen_set - set(metadata_by_run)),
            "duplicate_assignments": len(duplicated_runs),
            "orphan_runs": len(ena_set - frozen_set),
        },
        "sample_accounting": {
            "biological_samples": len(samples),
            "runs_per_sample_min": min((len(value) for value in run_groups.values()), default=0),
            "runs_per_sample_max": max((len(value) for value in run_groups.values()), default=0),
            "deseq2_observations": len(samples),
            "technical_runs_are_independent_observations": False,
            "aggregation_policy": AGGREGATION_POLICY,
        },
        "design": {
            "formula": formula,
            "variable": variable,
            "covariates": covariates,
            "interactions": [],
            "levels": levels,
            "reference_level": reference_level,
            "matrix_rows": len(matrix),
            "matrix_columns": len(design_columns),
            "rank": rank,
            "full_rank": full_rank,
            "contrasts": len(contrast_rows),
            "all_contrasts_estimable": all(row["estimable"] == "TRUE" for row in contrast_rows),
        },
        "ena_snapshot": {
            "source": "ENA Portal API",
            "query": next(iter(ena_queries), ""),
            "rows": len(ena),
            "sha256": sha256(ena_path),
        },
        "errors": errors,
    }
    (output / "preflight.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--ena-tsv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output_dir or root / "provenance" / EXPECTED_STUDY
    report = validate(root, args.ena_tsv, output)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
