"""Prepare species-specific BLAST checks for candidates absent from pooled top 100."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


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


def write(path: Path, records: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n")
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i:i + 60] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--forward", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with args.forward.open(encoding="utf-8", newline="") as handle:
        absent = [(row["gene_id"], row["species_code"]) for row in csv.DictReader(handle, delimiter="\t") if not row["best_hit_protein"]]
    query, panel = fasta(args.query), fasta(args.panel)
    for prefix in sorted({p for _, p in absent}):
        genes = {gene for gene, p in absent if p == prefix}
        write(args.output / f"{prefix}_query.faa", {name: sequence for name, sequence in query.items() if name.split(".")[0] in genes})
        write(args.output / f"{prefix}_reference.faa", {name: sequence for name, sequence in panel.items() if name.startswith(prefix + "|")})
    print(absent)


if __name__ == "__main__":
    main()
