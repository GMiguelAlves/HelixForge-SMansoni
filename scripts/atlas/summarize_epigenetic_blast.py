"""Summarize reproducible BLASTP results without treating hits as validation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "query_protein", "subject_protein", "identity_percent", "alignment_aa",
    "query_aa", "subject_aa", "query_start", "query_end", "subject_start",
    "subject_end", "evalue", "bit_score", "query_coverage_percent", "subject_title",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_hits(path: Path) -> dict[str, list[dict[str, str]]]:
    by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, fieldnames=FIELDS, delimiter="\t"):
            match = re.fullmatch(r"(Smp_\d+)\.\d+", row["query_protein"])
            if not match:
                raise ValueError(f"unexpected query protein: {row['query_protein']}")
            by_gene[match.group(1)].append(row)
    return by_gene


def span_coverage(intervals: list[tuple[int, int]], length: int) -> float:
    covered = set()
    for start, end in intervals:
        covered.update(range(min(start, end), max(start, end) + 1))
    return round(100 * len(covered) / length, 1)


def collapsed_hits(hits: list[dict[str, str]]) -> list[dict[str, object]]:
    """Rank subject/query pairs by combined HSP score; retain HSP identity separately."""
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for hit in hits:
        groups[(hit["query_protein"], hit["subject_protein"])].append(hit)
    result = []
    for group in groups.values():
        best_hsp = max(group, key=lambda row: float(row["bit_score"]))
        result.append({
            **best_hsp,
            "total_bit_score": round(sum(float(row["bit_score"]) for row in group), 1),
            "hsp_count": len(group),
            "query_coverage_percent": max(float(row["query_coverage_percent"]) for row in group),
            "subject_coverage_percent": span_coverage(
                [(int(row["subject_start"]), int(row["subject_end"])) for row in group],
                int(best_hsp["subject_aa"]),
            ),
        })
    return sorted(result, key=lambda row: float(row["total_bit_score"]), reverse=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query-manifest", type=Path, required=True)
    parser.add_argument("--human", type=Path, required=True)
    parser.add_argument("--self-blast", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.query_manifest.open(encoding="utf-8", newline="") as handle:
        queries = list(csv.DictReader(handle, delimiter="\t"))
    genes = list(dict.fromkeys(row["gene_id"] for row in queries))
    if len(genes) != 118:
        raise ValueError(f"expected 118 genes, found {len(genes)}")
    human = read_hits(args.human)
    self_hits = read_hits(args.self_blast)

    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for gene in genes:
        hits = collapsed_hits(human.get(gene, []))
        best = hits[0] if hits else None
        distinct_human = list(dict.fromkeys(row["subject_protein"] for row in hits))[:3]
        other = next(
            (row for row in collapsed_hits(self_hits.get(gene, [])) if not row["subject_protein"].startswith(gene + ".")),
            None,
        )
        other_gene = other["subject_protein"].split(".", 1)[0] if other else ""
        notes = []
        if not best:
            notes.append("NO_HUMAN_REVIEWED_HIT_AT_EVALUE_0.001")
        elif float(best["query_coverage_percent"]) < 50:
            notes.append("HUMAN_ALIGNMENT_COVERS_LESS_THAN_HALF_QUERY")
        if any(row["internal_stop"] == "YES" for row in queries if row["gene_id"] == gene):
            notes.append("INTERNAL_STOP_IN_PROTEIN")
        rows.append({
            "gene_id": gene,
            "isoforms": sum(row["gene_id"] == gene for row in queries),
            "best_query_isoform": best["query_protein"] if best else "",
            "best_human_accession": best["subject_protein"] if best else "",
            "best_human_title": best["subject_title"] if best else "",
            "best_human_best_hsp_identity_percent": best["identity_percent"] if best else "",
            "best_human_query_coverage_percent": best["query_coverage_percent"] if best else "",
            "best_human_subject_coverage_percent": best["subject_coverage_percent"] if best else "",
            "best_human_best_hsp_alignment_aa": best["alignment_aa"] if best else "",
            "best_human_hsp_count": best["hsp_count"] if best else "",
            "best_human_evalue": best["evalue"] if best else "",
            "best_human_total_bit_score": best["total_bit_score"] if best else "",
            "top_three_human_accessions": ";".join(distinct_human),
            "nearest_sm_paralog": other_gene,
            "nearest_sm_paralog_identity_percent": other["identity_percent"] if other else "",
            "nearest_sm_paralog_query_coverage_percent": other["query_coverage_percent"] if other else "",
            "review_flags": ";".join(notes),
        })

    output_path = args.output / "blastp_first_pass_by_gene.tsv"
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "interpretation": "sequence_similarity_only_not_functional_validation",
        "genes": len(rows),
        "protein_isoforms": len(queries),
        "genes_with_human_reviewed_hit_evalue_0_001": sum(bool(row["best_human_accession"]) for row in rows),
        "genes_with_human_alignment_under_50_percent_query_coverage": sum("HUMAN_ALIGNMENT_COVERS_LESS_THAN_HALF_QUERY" in row["review_flags"] for row in rows),
        "genes_with_nonself_sm_paralog_hit": sum(bool(row["nearest_sm_paralog"]) for row in rows),
        "sha256": {"query_manifest": sha256(args.query_manifest), "human_blast": sha256(args.human), "self_blast": sha256(args.self_blast), "summary_table": sha256(output_path)},
    }
    (args.output / "blastp_first_pass_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
