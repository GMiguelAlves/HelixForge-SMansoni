"""Map published corrected sirtuin proteins to current SM_V10 protein models."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from summarize_epigenetic_blast import FIELDS, collapsed_hits


def fasta(path: Path) -> dict[str, str]:
    records = {}
    name, parts = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name:
                records[name] = "".join(parts)
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    if name:
        records[name] = "".join(parts)
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--published", type=Path, required=True)
    parser.add_argument("--smansoni", type=Path, required=True)
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    published, sm = fasta(args.published), fasta(args.smansoni)
    hits = defaultdict(list)
    with args.blast.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, fieldnames=FIELDS, delimiter="\t"):
            hits[row["query_protein"]].append(row)
    data = []
    for accession, reference in published.items():
        ranked = collapsed_hits(hits.get(accession, []))
        best = ranked[0] if ranked else None
        current_id = best["subject_protein"] if best else ""
        current = sm.get(current_id, "")
        if current and len(current) <= len(reference):
            suffix = reference[-len(current):]
            suffix_mismatches = sum(a != b for a, b in zip(current, suffix))
        else:
            suffix_mismatches = None
        close_suffix = suffix_mismatches is not None and suffix_mismatches <= 5
        data.append({"published_accession": accession, "published_length_aa": len(reference),
            "best_sm_v10_protein": current_id, "best_sm_v10_gene": current_id.split(".")[0] if current_id else "",
            "sm_v10_length_aa": len(current) if current else "",
            "best_hsp_identity_percent": best["identity_percent"] if best else "",
            "published_query_coverage_percent": best["query_coverage_percent"] if best else "",
            "sm_v10_subject_coverage_percent": best["subject_coverage_percent"] if best else "",
            "evalue": best["evalue"] if best else "",
            "published_suffix_mismatches_if_at_most_5": suffix_mismatches if close_suffix else "",
            "missing_n_terminal_aa_if_suffix_match": len(reference) - len(current) if close_suffix else "",
        })
    path = args.output / "published_sirtuins_vs_sm_v10.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(data)
    (args.output / "published_sirtuins_mapping_manifest.json").write_text(json.dumps({
        "proteins": len(data), "published_sha256": hashlib.sha256(args.published.read_bytes()).hexdigest(),
        "sm_v10_sha256": hashlib.sha256(args.smansoni.read_bytes()).hexdigest(),
        "blast_sha256": hashlib.sha256(args.blast.read_bytes()).hexdigest(),
        "table_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(data, indent=2))


if __name__ == "__main__":
    main()
