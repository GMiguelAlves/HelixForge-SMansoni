"""Summarize top same-genus hits and prepare reciprocal BLAST queries."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from summarize_epigenetic_blast import collapsed_hits, read_hits


PREFIX = {"SHAE": "Schistosoma haematobium", "SJAP": "Schistosoma japonicum", "SROD": "Schistosoma rodhaini"}


def fasta(path: Path) -> dict[str, str]:
    result = {}
    name, pieces = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name:
                result[name] = "".join(pieces)
            name, pieces = line[1:].split()[0], []
        else:
            pieces.append(line.strip())
    if name:
        result[name] = "".join(pieces)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    genes = list(dict.fromkeys(x.split(".")[0] for x in fasta(args.queries)))
    hits = read_hits(args.blast)
    table = []
    best_subjects = set()
    for gene in genes:
        collapsed = collapsed_hits(hits.get(gene, []))
        for prefix, species in PREFIX.items():
            best = next((row for row in collapsed if row["subject_protein"].startswith(prefix + "|")), None)
            if best:
                best_subjects.add(str(best["subject_protein"]))
            table.append({
                "gene_id": gene, "species": species, "species_code": prefix,
                "query_isoform": best["query_protein"] if best else "",
                "best_hit_protein": best["subject_protein"] if best else "",
                "best_hsp_identity_percent": best["identity_percent"] if best else "",
                "query_coverage_percent": best["query_coverage_percent"] if best else "",
                "subject_coverage_percent": best["subject_coverage_percent"] if best else "",
                "evalue": best["evalue"] if best else "",
                "total_bit_score": best["total_bit_score"] if best else "",
                "hsp_count": best["hsp_count"] if best else "",
                "interpretation": "similarity_not_orthology" if best else "not_in_top_100_hits",
            })
    with (args.output / "best_close_species_hits.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(table)
    reference = fasta(args.reference)
    missing = best_subjects - set(reference)
    if missing:
        raise ValueError(f"missing {len(missing)} selected reference proteins: {list(missing)[:3]}")
    with (args.output / "best_close_species_queries.faa").open("w", encoding="utf-8") as handle:
        for subject in sorted(best_subjects):
            handle.write(f">{subject}\n")
            sequence = reference[subject]
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index:index + 60] + "\n")
    manifest = {"genes": len(genes), "rows": len(table), "selected_reference_isoforms": len(best_subjects),
                "genes_with_hit_by_species": {prefix: sum(r["species_code"] == prefix and bool(r["best_hit_protein"]) for r in table) for prefix in PREFIX},
                "blast_sha256": hashlib.sha256(args.blast.read_bytes()).hexdigest(),
                "reference_fasta_sha256": hashlib.sha256(args.reference.read_bytes()).hexdigest()}
    (args.output / "close_species_summary.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
