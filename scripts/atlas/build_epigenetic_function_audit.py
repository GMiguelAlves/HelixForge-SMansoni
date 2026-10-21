"""Combine first-pass epigenetic sequence evidence without inferring function."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
from collections import Counter, defaultdict
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, data: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(data)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--query-manifest", type=Path, required=True)
    parser.add_argument("--query-fasta", type=Path, required=True)
    parser.add_argument("--blast-summary", type=Path, required=True)
    parser.add_argument("--human-blast", type=Path, required=True)
    parser.add_argument("--self-blast", type=Path, required=True)
    parser.add_argument("--uniprot-domains", type=Path, required=True)
    parser.add_argument("--orthologs-gz", type=Path, required=True)
    parser.add_argument("--human-reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    selected = rows(args.selected)
    genes = list(dict.fromkeys(r["resolved_gene_id"] for r in selected if r["present_in_reference"] == "YES"))
    if len(genes) != 118:
        raise ValueError("expected exactly 118 distinct resolved genes")
    selected_by = {r["resolved_gene_id"]: r for r in selected if r["present_in_reference"] == "YES"}
    blast = {r["gene_id"]: r for r in rows(args.blast_summary)}
    domains = {r["gene_id"]: r for r in rows(args.uniprot_domains)}
    if set(genes) != set(blast) or set(genes) != set(domains):
        raise ValueError("gene IDs differ between selected, BLAST and domain inputs")

    same_genus: dict[str, list[dict[str, str]]] = defaultdict(list)
    ortholog_count = Counter()
    with gzip.open(args.orthologs_gz, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row["gene_id"]
            if gene not in selected_by:
                continue
            ortholog_count[gene] += 1
            if row["ortholog_species_name"].startswith("Schistosoma "):
                same_genus[gene].append(row)

    args.output.mkdir(parents=True, exist_ok=True)
    ortholog_data = [r for gene in genes for r in same_genus[gene]]
    if ortholog_data:
        write_tsv(args.output / "schistosoma_orthologs.tsv", ortholog_data)

    audit = []
    for gene in genes:
        source, hit, domain = selected_by[gene], blast[gene], domains[gene]
        same = same_genus[gene]
        flags = [f for f in hit["review_flags"].split(";") if f]
        if not domain["uniprot_length"]:
            flags.append("UNIPROT_RECORD_WITHOUT_LENGTH")
        if not domain["pfam"] and not domain["interpro"]:
            flags.append("NO_UNIPROT_PFAM_OR_INTERPRO_CROSS_REFERENCE")
        if not same:
            flags.append("NO_SCHISTOSOMA_ORTHOLOG_IN_WBPS19_FILE")
        audit.append({
            "gene_id": gene,
            "atlas_group": source["groups"],
            "atlas_functional_name": source["functional_name"],
            "protein_isoforms": hit["isoforms"],
            "best_blast_query_isoform": hit["best_query_isoform"],
            "best_human_reviewed_protein": hit["best_human_accession"],
            "best_human_title": hit["best_human_title"],
            "best_human_hsp_identity_percent": hit["best_human_best_hsp_identity_percent"],
            "best_human_query_coverage_percent": hit["best_human_query_coverage_percent"],
            "best_human_subject_coverage_percent": hit["best_human_subject_coverage_percent"],
            "best_human_evalue": hit["best_human_evalue"],
            "top_three_human_proteins": hit["top_three_human_accessions"],
            "nearest_sm_paralog": hit["nearest_sm_paralog"],
            "uniprot_accession": domain["source_accession"],
            "uniprot_protein_name": domain["protein_name"],
            "uniprot_length_aa": domain["uniprot_length"],
            "pfam_ids": domain["pfam"],
            "interpro_ids": domain["interpro"],
            "uniprot_feature_domains": domain["feature_domains"],
            "wbps19_ortholog_rows_all_species": str(ortholog_count[gene]),
            "wbps19_ortholog_rows_schistosoma": str(len(same)),
            "wbps19_schistosoma_species": ";".join(sorted({r["ortholog_species_name"] for r in same})),
            "review_flags": ";".join(flags),
            "interpretation_limit": "Sequence/domain evidence only; specific molecular function not experimentally validated by this audit",
        })
    write_tsv(args.output / "epigenetic_118_first_pass.tsv", audit)

    supporting = {
        "all_isoforms.faa": args.query_fasta,
        "isoform_manifest.tsv": args.query_manifest,
        "isoform_manifest.json": args.query_manifest.with_suffix(".json"),
        "human_blastp_raw.tsv": args.human_blast,
        "smansoni_blastp_raw.tsv": args.self_blast,
        "uniprot_domains.tsv": args.uniprot_domains,
    }
    for name, path in supporting.items():
        shutil.copyfile(path, args.output / name)

    summary = {
        "scope": "first_pass_sequence_evidence_not_functional_validation",
        "genes": len(genes),
        "protein_isoforms": sum(int(r["protein_isoforms"]) for r in audit),
        "genes_with_human_reviewed_blast_hit": sum(bool(r["best_human_reviewed_protein"]) for r in audit),
        "genes_with_best_human_query_coverage_below_50_percent": sum(
            bool(r["best_human_reviewed_protein"]) and float(r["best_human_query_coverage_percent"]) < 50
            for r in audit
        ),
        "genes_with_uniprot_pfam_or_interpro": sum(bool(r["pfam_ids"] or r["interpro_ids"]) for r in audit),
        "genes_with_schistosoma_ortholog_row": sum(bool(same_genus[g]) for g in genes),
        "sources": {
            "smansoni_proteins": "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS19.protein.fa.gz",
            "wbps19_orthologs": "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS19.orthologs.tsv.gz",
            "human_reviewed": "https://rest.uniprot.org/uniprotkb/stream?query=%28organism_id%3A9606%29%20AND%20%28reviewed%3Atrue%29&format=fasta",
            "uniprot_entries": "https://rest.uniprot.org/uniprotkb/search",
        },
        "blastp": {"version": "2.17.0+", "evalue_max": 0.001, "max_hsps": 10, "max_target_seqs_human": 20, "max_target_seqs_smansoni": 30},
        "input_sha256": {key: sha256(getattr(args, key)) for key in (
            "selected", "query_manifest", "query_fasta", "blast_summary", "human_blast",
            "self_blast", "uniprot_domains", "orthologs_gz", "human_reference",
        )},
        "output_sha256": {p.name: sha256(p) for p in sorted(args.output.iterdir()) if p.is_file() and p.suffix != ".json"},
    }
    (args.output / "run_manifest.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in summary.items() if key.startswith("genes") or key == "protein_isoforms"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
