"""Compare parasite Pfam models with each best reviewed-human BLAST hit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-pass", type=Path, required=True)
    parser.add_argument("--parasite-domains", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    first = read(args.first_pass)
    domains = defaultdict(set)
    for row in read(args.parasite_domains):
        domains[row["gene_id"]].update(x for x in row["pfam_ids"].split(";") if x)
    accession_by_gene = {}
    for row in first:
        match = re.match(r"sp\|([^|]+)\|", row["best_human_reviewed_protein"])
        if match:
            accession_by_gene[row["gene_id"]] = match.group(1)
    accessions = list(dict.fromkeys(accession_by_gene.values()))
    human = {}
    urls = []
    for index in range(0, len(accessions), 30):
        batch = accessions[index:index + 30]
        query = "(" + " OR ".join(f"accession:{accession}" for accession in batch) + ")"
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({
            "query": query, "format": "tsv", "fields": "accession,protein_name,length,xref_pfam", "size": 500})
        urls.append(url)
        for attempt in range(3):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "helixforge-epigenetic-audit/1.0"})
                payload = urllib.request.urlopen(request, timeout=60).read().decode("utf-8")
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(attempt + 1)
        (args.output / f"human_pfam_batch_{index // 30 + 1}.tsv").write_text(payload, encoding="utf-8")
        for row in csv.DictReader(io.StringIO(payload), delimiter="\t"):
            human[row["Entry"]] = row
    result = []
    for row in first:
        gene = row["gene_id"]
        accession = accession_by_gene.get(gene, "")
        reference = human.get(accession, {})
        reference_pfams = {p for p in reference.get("Pfam", "").split(";") if p}
        parasite_pfams = domains[gene]
        shared = sorted(parasite_pfams & reference_pfams)
        if not accession:
            status = "NO_HUMAN_BLAST_HIT"
        elif accession not in human:
            status = "HUMAN_ACCESSION_NOT_RETURNED"
        elif not parasite_pfams or not reference_pfams:
            status = "NO_PFAM_SET_ON_ONE_SIDE"
        elif shared:
            status = "SHARED_PFAM_MODEL"
        else:
            status = "NO_SHARED_PFAM_MODEL"
        result.append({"gene_id": gene, "best_human_accession": accession,
            "human_protein_name": reference.get("Protein names", ""),
            "human_pfam_ids": ";".join(sorted(reference_pfams)),
            "parasite_pfam_ids": ";".join(sorted(parasite_pfams)),
            "shared_pfam_ids": ";".join(shared),
            "comparison_status": status,
            "interpretation": "Pfam set overlap is family-level context, not substrate or orthology proof"})
    table = args.output / "human_parasite_pfam_comparison.tsv"
    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(result)
    from collections import Counter
    summary = {"genes": len(result), "status_counts": dict(Counter(r["comparison_status"] for r in result)),
               "source_queries": urls, "first_pass_sha256": sha(args.first_pass),
               "parasite_domains_sha256": sha(args.parasite_domains), "table_sha256": sha(table)}
    (args.output / "human_domain_comparison_manifest.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary["status_counts"], indent=2))


if __name__ == "__main__":
    main()
