"""Prepare a focused SIRT1-7/Schistosoma alignment panel for paralog review."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


def fasta(path: Path) -> dict[str, tuple[str, str]]:
    result = {}
    header, sequence = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                result[header.split()[0]] = (header, "".join(sequence).rstrip("*"))
            header, sequence = line[1:], []
        else:
            sequence.append(line.strip())
    if header:
        result[header.split()[0]] = (header, "".join(sequence).rstrip("*"))
    return result


def write_fasta(path: Path, data: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in data.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--human", type=Path, required=True)
    parser.add_argument("--smansoni", type=Path, required=True)
    parser.add_argument("--close-species", type=Path, required=True)
    parser.add_argument("--close-hits", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    human = fasta(args.human)
    human_sirtuins = {}
    for name, (header, sequence) in human.items():
        match = re.search(r"\bGN=(SIRT[1-7])\b", header)
        if match:
            human_sirtuins[match.group(1)] = (name, sequence)
    if set(human_sirtuins) != {f"SIRT{x}" for x in range(1, 8)}:
        raise ValueError(f"missing human sirtuins: {set(human_sirtuins)}")
    write_fasta(args.output / "human_sirtuins.faa", {f"HUMAN_{gene}": pair[1] for gene, pair in sorted(human_sirtuins.items())})

    # Include the five S. mansoni sirtuin-family candidates, including SIRT1/5
    # outside the atlas collection. A reciprocal hit to an unrelated gene is
    # not a reason to add that gene to this domain-specific tree.
    mansoni_genes = {"Smp_084140", "Smp_210340", "Smp_340470", "Smp_138640", "Smp_055090"}
    sm = fasta(args.smansoni)
    sm_by_gene = {}
    for protein, (_, sequence) in sm.items():
        gene = protein.split(".")[0]
        if gene in mansoni_genes and (gene not in sm_by_gene or len(sequence) > len(sm_by_gene[gene][1])):
            sm_by_gene[gene] = (protein, sequence)
    if set(sm_by_gene) != mansoni_genes:
        raise ValueError(f"missing SM_V10 sirtuin candidates: {mansoni_genes - set(sm_by_gene)}")
    close = fasta(args.close_species)
    with args.close_hits.open(encoding="utf-8", newline="") as handle:
        hits = [r for r in csv.DictReader(handle, delimiter="\t")
                if r["gene_id"] in {"Smp_084140", "Smp_210340", "Smp_340470"}
                and r["reciprocal_best_matches_original"] == "YES"]
    all_sequences = {f"HUMAN_{gene}": pair[1] for gene, pair in sorted(human_sirtuins.items())}
    all_sequences.update({gene: pair[1] for gene, pair in sorted(sm_by_gene.items())})
    for row in hits:
        protein = row["best_hit_protein"]
        if protein:
            name = row["species_code"] + "_" + row["gene_id"]
            all_sequences[name] = close[protein][1]
    write_fasta(args.output / "sirtuins_full_length.faa", all_sequences)
    manifest = {"human_uniprot": {gene: pair[0] for gene, pair in human_sirtuins.items()},
                "smansoni_proteins": {gene: pair[0] for gene, pair in sm_by_gene.items()},
                "close_hit_proteins": {r["species_code"] + "_" + r["gene_id"]: r["best_hit_protein"] for r in hits if r["best_hit_protein"]},
                "sequences": len(all_sequences), "sequence_sha256": hashlib.sha256((args.output / "sirtuins_full_length.faa").read_bytes()).hexdigest()}
    (args.output / "sirtuin_input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"sequences": len(all_sequences), "human": len(human_sirtuins), "smansoni": len(sm_by_gene), "close": len(hits)}, indent=2))


if __name__ == "__main__":
    main()
