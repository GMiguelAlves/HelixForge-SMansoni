"""Extract every WBPS19 protein isoform of the 118 atlas epigenetic genes."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_records(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        header = None
        sequence = []
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header, sequence = line[1:], []
            elif line:
                sequence.append(line)
        if header is not None:
            yield header, "".join(sequence)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--proteins", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.selected.open(encoding="utf-8", newline="") as handle:
        selected = list(csv.DictReader(handle, delimiter="\t"))
    genes = list(dict.fromkeys(row["resolved_gene_id"] for row in selected if row["present_in_reference"] == "YES"))
    if len(genes) != 118:
        raise ValueError(f"expected 118 unique atlas genes, found {len(genes)}")

    hits: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for header, sequence in fasta_records(args.proteins):
        match = re.search(r"(?:^|\s)gene=(Smp_\d+)(?:\s|$)", header)
        if match and match.group(1) in genes:
            hits[match.group(1)].append((header.split()[0], sequence))

    args.output.mkdir(parents=True, exist_ok=True)
    query_path = args.output / "epigenetic_118_all_isoforms.faa"
    manifest_path = args.output / "query_manifest.tsv"
    with query_path.open("w", encoding="utf-8") as fasta, manifest_path.open("w", encoding="utf-8", newline="") as table:
        writer = csv.writer(table, delimiter="\t")
        writer.writerow(["gene_id", "protein_id", "length_aa", "internal_stop", "sequence_sha256"])
        for gene in genes:
            for protein, sequence in hits[gene]:
                fasta.write(f">{protein} gene={gene}\n")
                for offset in range(0, len(sequence), 60):
                    fasta.write(sequence[offset:offset + 60] + "\n")
                writer.writerow([gene, protein, len(sequence), "YES" if "*" in sequence.rstrip("*") else "NO", hashlib.sha256(sequence.encode()).hexdigest()])

    missing = [gene for gene in genes if not hits[gene]]
    result = {
        "selected_genes_sha256": sha256(args.selected),
        "wbps19_proteins_sha256": sha256(args.proteins),
        "queries_sha256": sha256(query_path),
        "genes_expected": len(genes),
        "genes_with_protein": len(genes) - len(missing),
        "protein_isoforms": sum(len(v) for v in hits.values()),
        "missing_genes": missing,
    }
    (args.output / "query_manifest.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
