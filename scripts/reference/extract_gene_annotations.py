#!/usr/bin/env python3
"""Extract a compact, deterministic gene annotation table from WBPS GFF3."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from urllib.parse import unquote


FIELDS = (
    "gene_id",
    "display_label",
    "functional_name",
    "annotation_name",
    "description",
    "biotype",
    "description_source",
    "source_accession",
    "previous_stable_id",
    "chromosome",
    "start",
    "end",
    "strand",
)


def parse_attributes(raw: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for item in raw.rstrip(";").split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        values[key] = unquote(value)
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows: dict[str, dict[str, str]] = {}
    with args.gff.open(encoding="utf-8") as handle:
        for raw_line in handle:
            if raw_line.startswith("#"):
                continue
            columns = raw_line.rstrip("\n").split("\t")
            if len(columns) != 9 or columns[2] != "gene":
                continue
            attributes = parse_attributes(columns[8])
            gene_id = attributes.get("ID", "").removeprefix("gene:")
            if not gene_id:
                continue
            annotation_name = attributes.get("Name", "")
            description = attributes.get("description", "")
            functional_name = (
                annotation_name
                if annotation_name and annotation_name != gene_id
                else description
            )
            if not functional_name:
                functional_name = "No functional description in WBPS19"
            rows[gene_id] = {
                "gene_id": gene_id,
                "display_label": f"{gene_id} — {functional_name}",
                "functional_name": functional_name,
                "annotation_name": annotation_name,
                "description": description,
                "biotype": attributes.get("biotype", ""),
                "description_source": attributes.get("description_source", ""),
                "source_accession": attributes.get("description_source_acc", ""),
                "previous_stable_id": attributes.get("previous_stable_id", ""),
                "chromosome": columns[0],
                "start": columns[3],
                "end": columns[4],
                "strand": columns[6],
            }

    if not rows:
        raise SystemExit("No gene records were extracted from the GFF3")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows[gene_id] for gene_id in sorted(rows))
    print(f"genes={len(rows)} output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
