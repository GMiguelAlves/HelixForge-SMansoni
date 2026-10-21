"""Add reciprocal SM_V10 search evidence to same-genus candidate homologs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from summarize_epigenetic_blast import FIELDS, collapsed_hits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--forward", type=Path, required=True)
    parser.add_argument("--reciprocal", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with args.forward.open(encoding="utf-8", newline="") as handle:
        forward = list(csv.DictReader(handle, delimiter="\t"))
    by_query = defaultdict(list)
    with args.reciprocal.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, fieldnames=FIELDS, delimiter="\t"):
            by_query[row["query_protein"]].append(row)
    result = []
    for row in forward:
        protein = row["best_hit_protein"]
        hits = collapsed_hits(by_query.get(protein, [])) if protein else []
        genes = list(dict.fromkeys(hit["subject_protein"].split(".")[0] for hit in hits))
        target = row["gene_id"]
        rank = genes.index(target) + 1 if target in genes else ""
        result.append({**row,
            "reciprocal_best_sm_gene": genes[0] if genes else "",
            "reciprocal_original_gene_rank": rank,
            "reciprocal_best_matches_original": "YES" if genes and genes[0] == target else "NO" if genes else "UNAVAILABLE",
        })
    output_path = args.output / "close_species_reciprocal_by_gene.tsv"
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(result)
    summary = {"rows": len(result), "reciprocal_best_matches": sum(r["reciprocal_best_matches_original"] == "YES" for r in result),
               "rows_without_forward_hit": sum(not r["best_hit_protein"] for r in result),
               "by_species": {species: {"forward_hits": sum(r["species_code"] == species and bool(r["best_hit_protein"]) for r in result),
                                       "reciprocal_best_matches": sum(r["species_code"] == species and r["reciprocal_best_matches_original"] == "YES" for r in result)}
                              for species in sorted({r["species_code"] for r in result})},
               "reciprocal_sha256": hashlib.sha256(args.reciprocal.read_bytes()).hexdigest()}
    (args.output / "close_species_reciprocal_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
