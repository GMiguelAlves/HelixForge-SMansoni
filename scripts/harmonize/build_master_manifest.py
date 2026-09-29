#!/usr/bin/env python3
"""Build a conservative run-level index and controlled vocabularies."""

from __future__ import annotations

import csv
from pathlib import Path


STUDIES = ("PRJNA602528", "PRJNA597909", "PRJEB14695", "PRJEB32839")
MASTER_FIELDS = (
    "study", "sample", "run", "layout", "stage", "sex", "tissue",
    "condition", "treatment", "time", "replicate", "batch", "source",
)
VOCABULARIES = ("stage", "sex", "tissue", "treatment")


def read(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    master: list[dict[str, str]] = []
    for study in STUDIES:
        metadata_dir = root / "metadata" / study
        metadata = read(metadata_dir / "metadata.csv", ",")
        runs = {row["run_accession"]: row for row in read(metadata_dir / "runs.tsv", "\t")}
        for row in metadata:
            run = runs[row["run_accession"]]
            master.append({
                "study": study,
                "sample": row["sample_id"],
                "run": row["run_accession"],
                "layout": run["library_layout"],
                "stage": row["stage"],
                "sex": row["sex"],
                "tissue": row["tissue"],
                "condition": row["condition"],
                "treatment": row["treatment"],
                "time": row["time_hours"],
                "replicate": row["replicate"],
                "batch": row["batch"],
                "source": "ENA",
            })

    master.sort(key=lambda row: (STUDIES.index(row["study"]), row["sample"], row["run"]))
    write_tsv(root / "metadata" / "master_sample_manifest.tsv", MASTER_FIELDS, master)

    vocab_dir = root / "metadata" / "vocabularies"
    for field in VOCABULARIES:
        values = sorted({row[field] for row in master})
        write_tsv(
            vocab_dir / f"{field}_vocabulary.tsv",
            ("value", "status", "note"),
            [
                {
                    "value": value,
                    "status": "source_preserved",
                    "note": "Initial registry value; no cross-study biological equivalence asserted.",
                }
                for value in values
            ],
        )


if __name__ == "__main__":
    main()

