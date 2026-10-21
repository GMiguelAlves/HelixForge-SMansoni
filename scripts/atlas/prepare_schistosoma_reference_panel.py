"""Fetch a version-matched three-species Schistosoma protein BLAST panel."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import urllib.request
from pathlib import Path


BASE = "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/"
SOURCES = {
    "SHAE": ("schistosoma_haematobium", "PRJNA78265"),
    "SJAP": ("schistosoma_japonicum", "PRJNA520774"),
    "SROD": ("schistosoma_rodhaini", "TD1_PRJEB44434"),
}


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    combined = args.output / "schistosoma_3species_wbps19.faa"
    manifest = {"species": {}, "combined_fasta": combined.name}
    with combined.open("w", encoding="utf-8") as target:
        for prefix, (species, project) in SOURCES.items():
            name = f"{species}.{project}.WBPS19.protein.fa.gz"
            url = f"{BASE}{species}/{project}/{name}"
            path = args.output / name
            if not path.exists():
                request = urllib.request.Request(url, headers={"User-Agent": "helixforge-epigenetic-audit/1.0"})
                with urllib.request.urlopen(request, timeout=120) as source, path.open("wb") as sink:
                    while block := source.read(1024 * 1024):
                        sink.write(block)
            count = 0
            with gzip.open(path, "rt", encoding="utf-8") as source:
                for line in source:
                    if line.startswith(">"):
                        count += 1
                        target.write(f">{prefix}|{line[1:].strip()}\n")
                    else:
                        target.write(line)
            manifest["species"][prefix] = {
                "species": species, "project": project, "url": url,
                "file": name, "compressed_sha256": sha256(path), "protein_isoforms": count,
            }
            print(f"{species}: {count} isoforms", flush=True)
    manifest["combined_sha256"] = sha256(combined)
    (args.output / "reference_panel_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
