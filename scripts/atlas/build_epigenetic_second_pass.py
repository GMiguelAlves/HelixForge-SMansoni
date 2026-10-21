"""Join direct Pfam evidence and close-species reciprocal searches for 118 genes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, data: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(data)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-pass", type=Path, required=True)
    parser.add_argument("--isoforms", type=Path, required=True)
    parser.add_argument("--exact-domains", type=Path, required=True)
    parser.add_argument("--direct-domains", type=Path, required=True)
    parser.add_argument("--close-reciprocal", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    first = read(args.first_pass)
    isoforms = read(args.isoforms)
    exact = read(args.exact_domains)
    direct = read(args.direct_domains)
    close = read(args.close_reciprocal)
    genes = [r["gene_id"] for r in first]
    if len(genes) != 118 or len(set(genes)) != 118 or len(isoforms) != 130 or len(close) != 354:
        raise ValueError("unexpected gene, isoform or close-species row count")
    by_isoform = defaultdict(list)
    for row in exact:
        by_isoform[row["protein_id"]].append({"pfam_id": row["pfam_id"], "name": row["pfam_name"],
            "start": int(row["start_aa"]), "end": int(row["end_aa"]), "method": "INTERPRO_EXACT_SEQUENCE"})
    for row in direct:
        by_isoform[row["protein_id"]].append({"pfam_id": row["pfam_id"], "name": row["pfam_name"],
            "start": int(row["start_aa"]), "end": int(row["end_aa"]), "method": "HMMER_PFAM_SCAN"})
    isoform_data = []
    for row in isoforms:
        protein = row["protein_id"]
        domains = sorted(by_isoform[protein], key=lambda x: (x["start"], x["end"], x["pfam_id"]))
        isoform_data.append({"gene_id": row["gene_id"], "protein_id": protein, "length_aa": row["length_aa"],
            "sequence_sha256": row["sequence_sha256"], "pfam_match_rows": str(len(domains)),
            "pfam_ids": ";".join(sorted({d["pfam_id"] for d in domains})),
            "pfam_locations": ";".join(f"{d['pfam_id']}:{d['start']}-{d['end']}" for d in domains),
            "domain_evidence_method": ";".join(sorted({d["method"] for d in domains}))})
    write(args.output / "domain_architecture_by_isoform.tsv", isoform_data)

    by_gene_isoforms = defaultdict(list)
    by_gene_close = defaultdict(list)
    for row in isoform_data:
        by_gene_isoforms[row["gene_id"]].append(row)
    for row in close:
        by_gene_close[row["gene_id"]].append(row)
    result = []
    for row in first:
        gene = row["gene_id"]
        members, relatives = by_gene_isoforms[gene], by_gene_close[gene]
        if len(relatives) != 3:
            raise ValueError(f"expected three species rows for {gene}")
        pfams = sorted({p for member in members for p in member["pfam_ids"].split(";") if p})
        hit_count = sum(bool(r["best_hit_protein"]) for r in relatives)
        reciprocal_count = sum(r["reciprocal_best_matches_original"] == "YES" for r in relatives)
        alternatives = sorted({r["reciprocal_best_sm_gene"] for r in relatives
                               if r["reciprocal_best_sm_gene"] and r["reciprocal_best_sm_gene"] != gene})
        coverages = [float(r["query_coverage_percent"]) for r in relatives if r["best_hit_protein"]]
        flags = [f for f in row["review_flags"].split(";") if f]
        if not pfams:
            flags.append("NO_PFAM_HIT_IN_THIS_ANALYSIS")
        if alternatives:
            flags.append("CLOSE_HOMOLOG_RECIPROCAL_BEST_IS_OTHER_SM_GENE")
        if hit_count < 3:
            flags.append("CLOSE_HOMOLOG_ABSENT_FROM_TOP_100_OR_EVALUE_FILTER")
        if coverages and min(coverages) < 50:
            flags.append("ONE_CLOSE_SPECIES_ALIGNMENT_COVERS_LESS_THAN_HALF_QUERY")
        if any(not member["pfam_ids"] for member in members):
            flags.append("AT_LEAST_ONE_ISOFORM_WITHOUT_PFAM_HIT")
        result.append({
            "gene_id": gene, "atlas_group": row["atlas_group"], "atlas_functional_name": row["atlas_functional_name"],
            "isoforms": row["protein_isoforms"], "pfam_ids_all_isoforms": ";".join(pfams),
            "isoforms_with_pfam_hit": str(sum(bool(member["pfam_ids"]) for member in members)),
            "best_human_reviewed_protein": row["best_human_reviewed_protein"],
            "best_human_query_coverage_percent": row["best_human_query_coverage_percent"],
            "best_human_subject_coverage_percent": row["best_human_subject_coverage_percent"],
            "close_species_with_hit": str(hit_count), "close_species_with_reciprocal_best": str(reciprocal_count),
            "minimum_close_species_query_coverage_percent": str(min(coverages)) if coverages else "",
            "alternative_sm_genes_in_reciprocal_search": ";".join(alternatives),
            "review_flags": ";".join(dict.fromkeys(flags)),
            "specific_molecular_function_demonstrated_by_this_analysis": "NO",
        })
    write(args.output / "epigenetic_118_second_pass.tsv", result)
    summary = {
        "genes": len(result), "protein_isoforms": len(isoform_data),
        "genes_with_pfam_hit": sum(bool(r["pfam_ids_all_isoforms"]) for r in result),
        "isoforms_with_pfam_hit": sum(bool(r["pfam_ids"]) for r in isoform_data),
        "close_species_forward_hits": sum(bool(r["best_hit_protein"]) for r in close),
        "close_species_reciprocal_best_matches": sum(r["reciprocal_best_matches_original"] == "YES" for r in close),
        "genes_with_reciprocal_best_match_in_all_three_species": sum(r["close_species_with_reciprocal_best"] == "3" for r in result),
        "genes_with_alternative_reciprocal_sm_gene": sum(bool(r["alternative_sm_genes_in_reciprocal_search"]) for r in result),
        "genes_without_pfam_hit": [r["gene_id"] for r in result if not r["pfam_ids_all_isoforms"]],
        "input_sha256": {name: sha(getattr(args, name)) for name in ("first_pass", "isoforms", "exact_domains", "direct_domains", "close_reciprocal")},
        "output_sha256": {p.name: sha(p) for p in args.output.glob("*.tsv")},
    }
    (args.output / "second_pass_manifest.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in summary.items() if key not in {"input_sha256", "output_sha256"}}, indent=2))


if __name__ == "__main__":
    main()
