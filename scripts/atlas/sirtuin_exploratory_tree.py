"""Exploratory sirtuin NJ tree with bootstrap; not an orthology call."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--vendor-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=100)
    args = parser.parse_args()
    sys.path.insert(0, str(args.vendor_dir.resolve()))
    import Bio  # noqa: PLC0415
    from Bio import AlignIO, Phylo  # noqa: PLC0415
    from Bio.Phylo.Consensus import bootstrap_trees, get_support  # noqa: PLC0415
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor  # noqa: PLC0415

    args.output.mkdir(parents=True, exist_ok=True)
    alignment = AlignIO.read(str(args.alignment), "fasta")
    calculator = DistanceCalculator("blosum62", skip_letters=("-", "X", "?"))
    constructor = DistanceTreeConstructor(calculator, "nj")
    distance = calculator.get_distance(alignment)
    tree = constructor.nj(distance)
    random.seed(42)
    replicates = list(bootstrap_trees(alignment, args.bootstrap, constructor))
    supported = get_support(tree, replicates, len(replicates))
    Phylo.write(supported, str(args.output / "sirtuins_nj_bootstrap.newick"), "newick")
    with (args.output / "sirtuin_human_distances.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["smansoni_gene", "human_sirtuin", "blosum62_alignment_distance"])
        for gene in ["Smp_084140", "Smp_210340", "Smp_340470", "Smp_138640", "Smp_055090"]:
            for human in [f"HUMAN_SIRT{x}" for x in range(1, 8)]:
                writer.writerow([gene, human, round(distance[gene, human], 5)])
    summary = {"method": "FAMSA alignment; >=70% occupancy trim; Biopython neighbor joining with BLOSUM62 distance",
        "bootstrap_replicates": args.bootstrap, "seed": 42, "biopython_version": Bio.__version__,
        "sequences": len(alignment), "columns": alignment.get_alignment_length(),
        "limitations": "Exploratory topology, full-length input and occupancy trim; no formal orthology inference",
        "alignment_sha256": hashlib.sha256(args.alignment.read_bytes()).hexdigest()}
    (args.output / "tree_manifest.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
