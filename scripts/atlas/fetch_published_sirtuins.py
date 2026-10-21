"""Fetch the experimentally corrected 2013 S. mansoni sirtuin sequences."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import urllib.request
from pathlib import Path


URLS = {
    "kc": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=KC993850,KC993851,KC993852,KC993853,KC993854,KC993855,KC993856,KC993857&rettype=fasta_cds_aa&retmode=text",
    "abg": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=ABG78545&rettype=fasta&retmode=text",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    sources = {}
    for key, url in URLS.items():
        request = urllib.request.Request(url, headers={"User-Agent": "helixforge-epigenetic-audit/1.0"})
        text = urllib.request.urlopen(request, timeout=60).read().decode("utf-8")
        path = args.output / f"published_sirtuins_{key}.faa"
        path.write_text(text, encoding="utf-8")
        sources[key] = {"url": url, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        for item in text.split(">"):
            if not item.strip():
                continue
            header, *sequence = item.splitlines()
            protein = re.search(r"protein_id=([^\]]+)", header)
            accession = protein.group(1) if protein else header.split()[0]
            gene = re.search(r"gene=(Sirt\d+)", header, re.IGNORECASE)
            gene_name = gene.group(1) if gene else "Sirt1"
            records.append((accession, gene_name, "".join(sequence).rstrip("*")))
    combined = args.output / "published_sirtuins_2013.faa"
    with combined.open("w", encoding="utf-8") as handle:
        for accession, gene, sequence in records:
            handle.write(f">{accession} gene={gene}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")
    manifest = {"source_paper": "https://doi.org/10.1371/journal.pntd.0002428",
                "proteins": [{"accession": acc, "gene": gene, "length_aa": len(seq)} for acc, gene, seq in records],
                "sources": sources, "combined_sha256": hashlib.sha256(combined.read_bytes()).hexdigest()}
    (args.output / "published_sirtuins_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest["proteins"], indent=2))


if __name__ == "__main__":
    main()
