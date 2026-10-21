"""Map InterPro Pfam matches only when UniProt and atlas protein sequences agree exactly."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    result = {}
    name, parts = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name:
                result[name] = "".join(parts).rstrip("*")
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    if name:
        result[name] = "".join(parts).rstrip("*")
    return result


def fetch(url: str) -> bytes:
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "helixforge-epigenetic-audit/1.0", "Accept": "application/json"})
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            if error.code == 204:
                return b""
            if error.code not in (429, 500, 502, 503, 504) or attempt == 3:
                raise
        except (TimeoutError, urllib.error.URLError):
            if attempt == 3:
                raise
        time.sleep(2 * (attempt + 1))
    raise RuntimeError(url)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    queries = read_fasta(args.queries)
    with args.selected.open(encoding="utf-8", newline="") as handle:
        selected = list(csv.DictReader(handle, delimiter="\t"))
    accession_by_gene = {row["resolved_gene_id"]: row["source_accession"] for row in selected if row["present_in_reference"] == "YES"}
    if len(accession_by_gene) != 118:
        raise ValueError("expected 118 selected genes")

    # UniProt returns FASTA for active entries; merged accessions have no sequence here.
    uniprot = {}
    urls = []
    accessions = list(dict.fromkeys(accession_by_gene.values()))
    for index in range(0, len(accessions), 30):
        batch = accessions[index:index + 30]
        query = "(" + " OR ".join(f"accession:{a}" for a in batch) + ")"
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({"query": query, "format": "fasta", "size": 500})
        urls.append(url)
        payload = fetch(url).decode("utf-8")
        path = args.output / f"uniprot_batch_{index // 30 + 1}.faa"
        path.write_text(payload, encoding="utf-8")
        uniprot.update(read_fasta(path))
    uniprot = {key.split("|")[1]: value for key, value in uniprot.items() if "|" in key}
    exact = {}
    for gene, accession in accession_by_gene.items():
        reference = uniprot.get(accession)
        exact[gene] = [protein for protein, sequence in queries.items() if protein.startswith(gene + ".") and sequence == reference]

    def retrieve(item: tuple[str, str]) -> tuple[str, dict | None, str]:
        gene, accession = item
        if not exact[gene]:
            return gene, None, "NO_EXACT_UNIPROT_SEQUENCE"
        url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{accession}/?page_size=200"
        payload = fetch(url)
        return gene, json.loads(payload) if payload else {"count": 0, "results": []}, url

    results = {}
    with ThreadPoolExecutor(max_workers=4) as executor:
        pending = {executor.submit(retrieve, item): item for item in accession_by_gene.items()}
        for future in as_completed(pending):
            gene, data, source = future.result()
            results[gene] = (data, source)
    table = []
    status = []
    for gene, accession in accession_by_gene.items():
        data, url = results[gene]
        matched = exact[gene]
        status.append({"gene_id": gene, "uniprot_accession": accession, "exact_isoforms": ";".join(matched),
                       "pfam_match_count": str(data.get("count", 0)) if data else "", "status": "EXACT" if matched else "NO_EXACT_UNIPROT_SEQUENCE"})
        if not data:
            continue
        if data.get("next"):
            raise ValueError(f"unexpected pagination for {gene}: {data['next']}")
        for result in data.get("results", []):
            meta = result["metadata"]
            for protein in result.get("proteins", []):
                for location in protein.get("entry_protein_locations", []):
                    for fragment in location.get("fragments", []):
                        for isoform in matched:
                            table.append({
                                "gene_id": gene, "protein_id": isoform, "uniprot_accession": accession,
                                "pfam_id": meta["accession"], "pfam_name": meta.get("name", ""),
                                "integrated_interpro": meta.get("integrated") or "",
                                "start_aa": fragment["start"], "end_aa": fragment["end"],
                                "score": location.get("score", ""), "source_url": url,
                            })
    for name, data in (("exact_sequence_status.tsv", status), ("exact_pfam_locations.tsv", table)):
        with (args.output / name).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(data)
    summary = {
        "genes": len(accession_by_gene),
        "genes_with_exact_uniprot_isoform": sum(bool(v) for v in exact.values()),
        "exact_isoforms": sum(map(len, exact.values())),
        "genes_with_pfam_locations": len({r["gene_id"] for r in table}),
        "source": "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{accession}/",
        "uniprot_fasta_queries": urls,
        "query_fasta_sha256": hashlib.sha256(args.queries.read_bytes()).hexdigest(),
    }
    (args.output / "exact_interpro_manifest.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in summary.items() if key.startswith("genes") or key == "exact_isoforms"}, indent=2))


if __name__ == "__main__":
    main()
