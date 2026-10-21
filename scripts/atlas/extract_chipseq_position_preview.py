#!/usr/bin/env python3
"""Extract a sanitized candidate-gene ChIP position table from an audit archive.

Only gene, stage, mark, normalized peak-position class, and peak genomic
coordinates are retained. Runtime paths and historical RNA-derived fields are
discarded. The exact-coordinate export is limited to the selected candidate
genes and is intended for the local atlas gene explorer.
"""

from __future__ import annotations

import argparse
import csv
import io
import tarfile
from collections import defaultdict
from pathlib import Path


def normalize_position(location: str, promoter: str) -> str:
    value = location.strip().lower().replace("-", "_").replace(" ", "_")
    if promoter.strip().lower() in {"true", "1", "yes"} or "promoter" in value or "tss" in value:
        return "promoter_TSS"
    if any(token in value for token in ("exon", "intron", "utr", "cds", "gene_body")):
        return "gene_body"
    if any(token in value for token in ("intergenic", "distal")):
        return "distal_intergenic"
    return "other_annotated"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--selected-genes", type=Path, required=True)
    parser.add_argument("--mark-stage-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--locus-output", type=Path)
    args = parser.parse_args()

    with args.selected_genes.open(encoding="utf-8", newline="") as handle:
        selected = {
            row["resolved_gene_id"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("resolved_gene_id")
        }
    with args.mark_stage_metadata.open(encoding="utf-8", newline="") as handle:
        valid_contexts = {
            (row["mark_or_factor"], row["stage_or_condition"])
            for row in csv.DictReader(handle, delimiter="\t")
        }

    grouped: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    locus_rows: dict[tuple[str, str, str, str], dict[str, str]] = {}
    member_name = "debug_integrateseq_tsv/gene_mark_stage_links.tsv"
    with tarfile.open(args.archive, "r:gz") as archive:
        member = archive.getmember(member_name)
        raw = archive.extractfile(member)
        if raw is None:
            raise RuntimeError(f"unable to read {member_name}")
        with io.TextIOWrapper(raw, encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                gene = row.get("gene_id", "")
                context = (row.get("mark_or_factor", ""), row.get("stage_or_condition", ""))
                if gene not in selected or context not in valid_contexts:
                    continue
                position = normalize_position(row.get("peak_location", ""), row.get("promoter_flag", ""))
                peak_id = row.get("peak_id", "")
                if peak_id:
                    grouped[(gene, context[1], context[0], position)].add(peak_id)
                    locus_rows[(gene, context[1], context[0], peak_id)] = {
                        "gene_id": gene,
                        "stage_or_condition": context[1],
                        "mark_or_factor": context[0],
                        "peak_id": peak_id,
                        "peak_chrom": row.get("peak_chrom", ""),
                        "peak_start": row.get("peak_start", ""),
                        "peak_end": row.get("peak_end", ""),
                        "peak_midpoint": row.get("peak_midpoint", ""),
                        "peak_location": row.get("peak_location", ""),
                        "promoter_flag": row.get("promoter_flag", ""),
                        "position_class": position,
                    }

    fields = ["gene_id", "stage_or_condition", "mark_or_factor", "position_class", "n_unique_peaks"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key in sorted(grouped):
            gene, stage, mark, position = key
            writer.writerow({
                "gene_id": gene,
                "stage_or_condition": stage,
                "mark_or_factor": mark,
                "position_class": position,
                "n_unique_peaks": len(grouped[key]),
            })
    print(f"wrote {len(grouped)} candidate gene/stage/mark/position records to {args.output}")
    if args.locus_output:
        locus_fields = [
            "gene_id", "stage_or_condition", "mark_or_factor", "peak_id",
            "peak_chrom", "peak_start", "peak_end", "peak_midpoint",
            "peak_location", "promoter_flag", "position_class",
        ]
        args.locus_output.parent.mkdir(parents=True, exist_ok=True)
        with args.locus_output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=locus_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for key in sorted(locus_rows):
                writer.writerow(locus_rows[key])
        print(f"wrote {len(locus_rows)} sanitized candidate peak loci to {args.locus_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
