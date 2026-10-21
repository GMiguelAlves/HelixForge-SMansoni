"""Index Europe PMC mentions of current SM_V10 IDs; not functional curation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path


BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def fetch(url: str) -> dict:
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "helixforge-epigenetic-audit/1.0"})
            with urllib.request.urlopen(request, timeout=45) as response:
                return json.load(response)
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(url)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with args.selected.open(encoding="utf-8", newline="") as handle:
        genes = list(dict.fromkeys(row["resolved_gene_id"] for row in csv.DictReader(handle, delimiter="\t")
                                   if row["present_in_reference"] == "YES"))
    if len(genes) != 118:
        raise ValueError(f"expected 118 genes; got {len(genes)}")
    summary, articles = [], []
    for index, gene in enumerate(genes, 1):
        url = BASE + "?" + urllib.parse.urlencode({"query": f'"{gene}"', "format": "json", "pageSize": 25, "resultType": "lite"})
        data = fetch(url)
        found = data.get("resultList", {}).get("result", [])
        summary.append({"gene_id": gene, "europe_pmc_hit_count": data.get("hitCount", 0),
                        "retrieved_first_n": len(found), "query_url": url,
                        "meaning": "ID mention only; no claim of direct functional evidence"})
        for item in found:
            articles.append({"gene_id": gene, "id": item.get("id", ""), "source": item.get("source", ""),
                             "pmcid": item.get("pmcid", ""), "doi": item.get("doi", ""),
                             "title": item.get("title", ""), "pub_year": item.get("pubYear", ""),
                             "query_url": url})
        if index % 20 == 0 or index == len(genes):
            print(f"searched {index}/{len(genes)} IDs", flush=True)
        time.sleep(0.5)
    for name, data in (("literature_mentions_by_gene.tsv", summary), ("literature_mention_articles.tsv", articles)):
        with (args.output / name).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(data)
    manifest = {"genes": len(genes), "genes_with_id_mentions": sum(int(r["europe_pmc_hit_count"]) > 0 for r in summary),
                "article_rows": len(articles), "source": BASE, "query_scope": "quoted current gene ID, Europe PMC search index",
                "limits": "IDs can be absent from abstracts/index despite being in supplement; mentions are not experimental validation",
                "selected_sha256": hashlib.sha256(args.selected.read_bytes()).hexdigest()}
    (args.output / "literature_search_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
