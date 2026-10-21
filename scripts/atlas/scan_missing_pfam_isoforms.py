"""Scan isoforms absent from exact UniProt matches against Pfam via EBI HMMER."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
import urllib.request
from pathlib import Path


BASE = "https://www.ebi.ac.uk/Tools/hmmer/api/v1"


def read_fasta(path: Path) -> dict[str, str]:
    records = {}
    name, chunks = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name:
                records[name] = "".join(chunks).rstrip("*")
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if name:
        records[name] = "".join(chunks).rstrip("*")
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--exact-status", type=Path, required=True)
    parser.add_argument("--exact-locations", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    queries = read_fasta(args.queries)
    with args.exact_status.open(encoding="utf-8", newline="") as handle:
        status = list(csv.DictReader(handle, delimiter="\t"))
    exact = {p for row in status for p in row["exact_isoforms"].split(";") if p}
    with args.exact_locations.open(encoding="utf-8", newline="") as handle:
        with_match = {row["protein_id"] for row in csv.DictReader(handle, delimiter="\t")}
    missing = {p: s for p, s in queries.items() if p not in exact or p not in with_match}
    table = []
    jobs = []
    for index, (protein, sequence) in enumerate(missing.items(), 1):
        result_path = args.output / f"{protein}.json"
        if result_path.exists():
            response = json.loads(result_path.read_text(encoding="utf-8"))
            job_id = response.get("result", {}).get("stats", {}).get("id", "")
        else:
            payload = json.dumps({"input": sequence, "database": "pfam", "threshold": "evalue", "E": 0.01,
                                  "domE": 0.01, "incE": 0.01, "incdomE": 0.01}).encode()
            request = urllib.request.Request(f"{BASE}/search/hmmscan", data=payload,
                headers={"Content-Type": "application/json", "User-Agent": "helixforge-epigenetic-audit/1.0"})
            with urllib.request.urlopen(request, timeout=90) as remote:
                job_id = json.load(remote)["id"]
            for attempt in range(50):
                with urllib.request.urlopen(f"{BASE}/result/{job_id}", timeout=60) as remote:
                    response = json.load(remote)
                if response.get("status") == "SUCCESS":
                    break
                if response.get("status") in {"FAILURE", "ERROR"}:
                    raise RuntimeError(f"{protein}: {response}")
                time.sleep(2)
            else:
                raise TimeoutError(f"{protein}: {job_id}")
            result_path.write_text(json.dumps(response, separators=(",", ":")) + "\n", encoding="utf-8")
            time.sleep(0.5)
        hits = response.get("result", {}).get("hits", [])
        for hit in hits:
            if not hit.get("is_included"):
                continue
            meta = hit.get("metadata") or {}
            for domain in hit.get("domains", []):
                if not domain.get("is_included") or domain.get("ievalue", 1) > 0.01:
                    continue
                table.append({"gene_id": protein.split(".")[0], "protein_id": protein,
                    "pfam_id": meta.get("accession", ""), "pfam_name": meta.get("identifier", ""),
                    "start_aa": domain["iali"], "end_aa": domain["jali"],
                    "domain_independent_evalue": domain["ievalue"], "domain_bitscore": domain["bitscore"],
                    "job_id": job_id})
        jobs.append({"protein_id": protein, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                     "job_id": job_id, "included_domain_rows": sum(r["protein_id"] == protein for r in table)})
        print(f"{index}/{len(missing)} {protein}: {jobs[-1]['included_domain_rows']} domain rows", flush=True)
    for name, data in (("hmmer_pfam_locations.tsv", table), ("hmmer_jobs.tsv", jobs)):
        with (args.output / name).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(data)
    (args.output / "scan_manifest.json").write_text(json.dumps({
        "method": "EBI HMMER hmmscan against Pfam, sequence and domain independent E <= 0.01",
        "endpoint": f"{BASE}/search/hmmscan", "queried_isoforms": len(missing),
        "isoforms_with_included_domain": len({r["protein_id"] for r in table}),
        "query_fasta_sha256": hashlib.sha256(args.queries.read_bytes()).hexdigest(),
    }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
