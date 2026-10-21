"""Fetch UniProt domain cross-references for the 118 SM_V10 atlas genes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--annotations", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    genes = list(dict.fromkeys(row["resolved_gene_id"] for row in read_tsv(args.selected) if row["present_in_reference"] == "YES"))
    if len(genes) != 118:
        raise ValueError(f"expected 118 unique genes, found {len(genes)}")
    annotations = {row["gene_id"]: row for row in read_tsv(args.annotations)}
    accessions = [annotations[gene]["source_accession"] for gene in genes]
    if any(not accession for accession in accessions) or len(set(accessions)) != 118:
        raise ValueError("missing or duplicate UniProt accessions")

    args.output.mkdir(parents=True, exist_ok=True)
    raw_rows: dict[str, dict[str, str]] = {}
    urls = []
    for offset in range(0, len(accessions), 30):
        batch = accessions[offset:offset + 30]
        query = "(" + " OR ".join(f"accession:{accession}" for accession in batch) + ")"
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({
            "query": query,
            "format": "tsv",
            "fields": "accession,id,protein_name,length,xref_pfam,xref_interpro,ft_domain",
            "size": 500,
        })
        urls.append(url)
        for attempt in range(3):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "helixforge-epigenetic-audit/1.0"})
                with urllib.request.urlopen(request, timeout=90) as response:
                    content = response.read().decode("utf-8")
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(2 * (attempt + 1))
        (args.output / f"uniprot_batch_{offset // 30 + 1}.tsv").write_text(content, encoding="utf-8")
        for row in csv.DictReader(io.StringIO(content), delimiter="\t"):
            raw_rows[row["Entry"]] = row

    columns = ["gene_id", "source_accession", "uniprot_entry", "protein_name", "uniprot_length", "pfam", "interpro", "feature_domains", "retrieval_status"]
    output_path = args.output / "uniprot_domains_by_gene.tsv"
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for gene, accession in zip(genes, accessions):
            row = raw_rows.get(accession)
            writer.writerow({
                "gene_id": gene,
                "source_accession": accession,
                "uniprot_entry": row.get("Entry Name", "") if row else "",
                "protein_name": row.get("Protein names", "") if row else "",
                "uniprot_length": row.get("Length", "") if row else "",
                "pfam": row.get("Pfam", "") if row else "",
                "interpro": row.get("InterPro", "") if row else "",
                "feature_domains": row.get("Domain [FT]", "") if row else "",
                "retrieval_status": "FOUND" if row else "NOT_FOUND",
            })

    summary = {
        "genes": len(genes),
        "uniprot_entries_found": sum(accession in raw_rows for accession in accessions),
        "uniprot_entries_missing": [gene for gene, accession in zip(genes, accessions) if accession not in raw_rows],
        "queries": urls,
        "sha256": {"selected": sha256(args.selected), "annotations": sha256(args.annotations), "domains_table": sha256(output_path)},
    }
    (args.output / "uniprot_domain_manifest.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: summary[key] for key in ("genes", "uniprot_entries_found", "uniprot_entries_missing")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
