"""Check the published 118-gene audit against the selected atlas collection."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    selected = {row["resolved_gene_id"] for row in rows(args.selected)}
    excluded = {row["sm_v10_gene_id"] for row in rows(args.candidates)
                if row["curation_status"].startswith("EXCLUDED_") and row["sm_v10_gene_id"]}
    first = rows(args.audit / "epigenetic_118_first_pass.tsv")
    second = rows(args.audit / "epigenetic_118_second_pass.tsv")
    isoforms = rows(args.audit / "isoform_manifest.tsv")
    architecture = rows(args.audit / "domain_architecture_by_isoform.tsv")
    close = rows(args.audit / "close_species_reciprocal_by_gene.tsv")
    human = rows(args.audit / "human_parasite_pfam_comparison.tsv")
    assert len(selected) == len(first) == len(second) == len(human) == 118
    assert selected == {row["gene_id"] for row in first} == {row["gene_id"] for row in second}
    assert selected == {row["gene_id"] for row in human}
    assert not selected & excluded
    assert len(isoforms) == len(architecture) == 130
    assert {row["protein_id"] for row in isoforms} == {row["protein_id"] for row in architecture}
    assert {row["gene_id"] for row in isoforms} == selected
    assert len(close) == 354
    assert set(Counter(row["gene_id"] for row in close).values()) == {3}
    assert Counter(row["comparison_status"] for row in human) == {
        "SHARED_PFAM_MODEL": 113, "NO_SHARED_PFAM_MODEL": 1,
        "NO_HUMAN_BLAST_HIT": 3, "NO_PFAM_SET_ON_ONE_SIDE": 1}
    published = {row["published_accession"]: row for row in rows(
        args.audit / "sirtuin_phylogeny" / "published_sirtuins_vs_sm_v10.tsv")}
    assert published["AGT95751.1"]["best_sm_v10_protein"] == "Smp_210340.1"
    assert published["AGT95751.1"]["published_suffix_mismatches_if_at_most_5"] == "0"
    assert published["AGT95752.1"]["best_sm_v10_protein"] == "Smp_340470.1"
    assert published["AGT95752.1"]["missing_n_terminal_aa_if_suffix_match"] == "119"
    assert published["AGT95752.1"]["published_suffix_mismatches_if_at_most_5"] == "1"
    manifest = json.loads((args.audit / "delivery_manifest.json").read_text(encoding="utf-8"))
    for name, expected in manifest["published_file_sha256"].items():
        actual = hashlib.sha256((args.audit / name).read_bytes()).hexdigest()
        assert actual == expected, f"hash mismatch: {name}"
    print(f"validated {len(selected)} genes, {len(isoforms)} isoforms, "
          f"{len(close)} species comparisons and {len(manifest['published_file_sha256'])} files")


if __name__ == "__main__":
    main()
