"""Publish reproducible second-pass inputs/results outside rebuilt atlas outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staging", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    files = {
        "interpro_exact/exact_sequence_status.tsv": "exact_sequence_status.tsv",
        "interpro_exact/exact_pfam_locations.tsv": "exact_pfam_locations.tsv",
        "interpro_exact/exact_interpro_manifest.json": "exact_interpro_manifest.json",
        "interpro_exact/uniprot_batch_1.faa": "exact_uniprot_reference_batch_1.faa",
        "interpro_exact/uniprot_batch_2.faa": "exact_uniprot_reference_batch_2.faa",
        "interpro_exact/uniprot_batch_3.faa": "exact_uniprot_reference_batch_3.faa",
        "interpro_exact/uniprot_batch_4.faa": "exact_uniprot_reference_batch_4.faa",
        "pfam_missing/hmmer_pfam_locations.tsv": "hmmer_pfam_locations.tsv",
        "pfam_missing/hmmer_jobs.tsv": "hmmer_jobs.tsv",
        "pfam_missing/scan_manifest.json": "hmmer_scan_manifest.json",
        "close_species/reference_panel_manifest.json": "close_species_reference_manifest.json",
        "close_species/close_species_blastp.tsv": "close_species_blastp_raw.tsv",
        "close_species/summary/best_close_species_hits.tsv": "best_close_species_hits.tsv",
        "close_species/summary/close_species_reciprocal_by_gene.tsv": "close_species_reciprocal_by_gene.tsv",
        "close_species/summary/close_species_reciprocal_summary.json": "close_species_reciprocal_summary.json",
        "data/close_species_reciprocal_blastp.tsv": "close_species_reciprocal_blastp_raw.tsv",
        "close_species/rescue/SHAE_rescue.tsv": "close_species_SHAE_rescue_blastp.tsv",
        "close_species/rescue/SJAP_rescue.tsv": "close_species_SJAP_rescue_blastp.tsv",
        "sirtuin_tree/sirtuins_full_length.faa": "sirtuin_phylogeny/sirtuins_full_length.faa",
        "sirtuin_tree/sirtuins_aligned.afa": "sirtuin_phylogeny/sirtuins_aligned.afa",
        "sirtuin_tree/sirtuins_70pct_occupancy.afa": "sirtuin_phylogeny/sirtuins_70pct_occupancy.afa",
        "sirtuin_tree/sirtuins_nj_bootstrap.newick": "sirtuin_phylogeny/sirtuins_nj_bootstrap.newick",
        "sirtuin_tree/sirtuin_human_distances.tsv": "sirtuin_phylogeny/sirtuin_human_distances.tsv",
        "sirtuin_tree/sirtuin_input_manifest.json": "sirtuin_phylogeny/sirtuin_input_manifest.json",
        "sirtuin_tree/alignment_manifest.json": "sirtuin_phylogeny/alignment_manifest.json",
        "sirtuin_tree/tree_manifest.json": "sirtuin_phylogeny/tree_manifest.json",
        "sirtuin_tree/published/published_sirtuins_2013.faa": "sirtuin_phylogeny/published_sirtuins_2013.faa",
        "sirtuin_tree/published/published_sirtuins_manifest.json": "sirtuin_phylogeny/published_sirtuins_manifest.json",
        "sirtuin_tree/published/published_sirtuins_vs_sm_v10.tsv": "sirtuin_phylogeny/published_sirtuins_vs_sm_v10.tsv",
        "sirtuin_tree/published/published_sirtuins_mapping_manifest.json": "sirtuin_phylogeny/published_sirtuins_mapping_manifest.json",
        "data/published_sirtuins_vs_smansoni.tsv": "sirtuin_phylogeny/published_sirtuins_vs_smansoni_blastp_raw.tsv",
        "literature/literature_mentions_by_gene.tsv": "literature_mentions_by_gene.tsv",
        "literature/literature_mention_articles.tsv": "literature_mention_articles.tsv",
        "literature/literature_search_manifest.json": "literature_search_manifest.json",
        "human_domain_comparison/human_parasite_pfam_comparison.tsv": "human_parasite_pfam_comparison.tsv",
        "human_domain_comparison/human_domain_comparison_manifest.json": "human_domain_comparison_manifest.json",
    }
    for batch in range(1, 5):
        files[f"human_domain_comparison/human_pfam_batch_{batch}.tsv"] = f"human_pfam_batch_{batch}.tsv"
    for source_name, target_name in files.items():
        source, target = args.staging / source_name, args.output / target_name
        if not source.exists():
            raise FileNotFoundError(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    raw_dir = args.output / "hmmer_raw"
    raw_dir.mkdir(exist_ok=True)
    with (args.staging / "pfam_missing/hmmer_jobs.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            source = args.staging / "pfam_missing" / f"{row['protein_id']}.json"
            shutil.copyfile(source, raw_dir / source.name)
    manifest = {"purpose": "second-pass computation, not automatic functional validation",
                "published_file_sha256": {str(path.relative_to(args.output)).replace("\\", "/"): sha256(path)
                     for path in sorted(args.output.rglob("*")) if path.is_file() and path.name != "delivery_manifest.json"},
                "upstream_sources": {
                    "WBPS19": "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/",
                    "InterPro": "https://www.ebi.ac.uk/interpro/api/",
                    "HMMER": "https://www.ebi.ac.uk/Tools/hmmer/api/v1/",
                    "published_sirtuins": "https://doi.org/10.1371/journal.pntd.0002428",
                }}
    (args.output / "delivery_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"published {len(manifest['published_file_sha256'])} files")


if __name__ == "__main__":
    main()
