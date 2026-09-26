#!/usr/bin/env python3
"""Validate the committed RNA-seq atlas bundle without optional dependencies."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


EXPECTED = {
    "studies": ["PRJNA602528", "PRJNA597909", "PRJEB14695", "PRJEB32839"],
    "biological_samples": 128,
    "technical_runs": 318,
    "genes": 9914,
    "contrasts": 22,
}
FORBIDDEN = (
    "/ho" + "me/",
    "/scr" + "atch/",
    "srv-" + "slurm",
    "@bio." + "ib.unicamp.br",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(root: Path) -> list[str]:
    errors: list[str] = []
    atlas = root / "results/atlas/rnaseq"
    manifest_path = atlas / "manifest.json"
    if not manifest_path.is_file():
        return ["atlas manifest is missing"]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for field, expected in EXPECTED.items():
        if manifest.get(field) != expected:
            errors.append(f"manifest {field}: expected {expected!r}, observed {manifest.get(field)!r}")
    if manifest.get("status") != "PASS":
        errors.append("atlas status is not PASS")
    if manifest.get("batch_effect_assessment") != "DEFERRED":
        errors.append("batch assessment must remain DEFERRED in atlas v1")
    if manifest.get("global_pca_interpretation") != "EXPLORATORY_NOT_A_BATCH_TEST":
        errors.append("global PCA interpretation boundary is missing")

    for record in manifest.get("inputs", []):
        path = root / record["path"]
        if not path.is_file():
            errors.append(f"missing input: {record['path']}")
            continue
        if sha256(path) != record["sha256"]:
            errors.append(f"input checksum mismatch: {record['path']}")

    for record in manifest.get("outputs", []):
        path = atlas / record["path"]
        if not path.is_file():
            errors.append(f"missing output: {record['path']}")
            continue
        if path.stat().st_size != record["bytes"]:
            errors.append(f"size mismatch: {record['path']}")
        if sha256(path) != record["sha256"]:
            errors.append(f"checksum mismatch: {record['path']}")

    samples = read_tsv(atlas / "data/sample_inventory.tsv")
    if len(samples) != EXPECTED["biological_samples"]:
        errors.append("sample inventory cardinality mismatch")
    if sum(int(row["technical_runs"]) for row in samples) != EXPECTED["technical_runs"]:
        errors.append("technical-run total mismatch")
    if {row["reference_id"] for row in samples} != {"Schistosoma_mansoni_SM_V10_WBPS19"}:
        errors.append("sample inventory mixes reference identities")
    if {row["batch"] for row in samples} != {"not_declared"}:
        errors.append("atlas v1 batch inventory no longer matches frozen metadata")

    contrasts = read_tsv(atlas / "data/contrast_catalog.tsv")
    if len(contrasts) != EXPECTED["contrasts"]:
        errors.append("contrast catalog cardinality mismatch")
    if {row["study"] for row in contrasts} != {"PRJNA597909", "PRJEB14695", "PRJEB32839"}:
        errors.append("unexpected studies in differential-expression catalog")

    for filename in ("gene_contrast_log2fc.tsv", "gene_contrast_padj.tsv"):
        path = atlas / "data" / filename
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            row_count = sum(1 for _ in reader)
        if len(header) - 1 != EXPECTED["contrasts"] or row_count != EXPECTED["genes"]:
            errors.append(f"unexpected dimensions: {filename}")

    payload_path = atlas / "data/atlas_payload.js"
    payload_text = payload_path.read_text(encoding="utf-8")
    if not payload_text.startswith("window.HF_ATLAS=") or not payload_text.endswith(";\n"):
        errors.append("browser payload wrapper is invalid")
    else:
        payload = json.loads(payload_text[len("window.HF_ATLAS="):-2])
        if len(payload.get("genes", [])) != EXPECTED["genes"]:
            errors.append("browser payload gene count mismatch")
        if len(payload.get("samples", [])) != EXPECTED["biological_samples"]:
            errors.append("browser payload sample count mismatch")
        if len(payload.get("contrasts", [])) != EXPECTED["contrasts"]:
            errors.append("browser payload contrast count mismatch")

    for svg in sorted((atlas / "figures").glob("*.svg")):
        try:
            root_element = ET.parse(svg).getroot()
            if not root_element.tag.endswith("svg"):
                errors.append(f"not an SVG document: {svg.name}")
        except ET.ParseError as exc:
            errors.append(f"invalid SVG {svg.name}: {exc}")

    report = (atlas / "atlas.html").read_text(encoding="utf-8")
    required_sections = ('id="overview"', 'id="samples"', 'id="biology"', 'id="genes"', 'id="data"')
    for marker in required_sections:
        if marker not in report:
            errors.append(f"missing report section: {marker}")
    for link in re.findall(r'(?:src|href)="([^"#]+)"', report):
        if "://" in link or link.startswith("mailto:"):
            continue
        if not (atlas / link).is_file():
            errors.append(f"broken local report link: {link}")

    for path in atlas.rglob("*"):
        if not path.is_file() or path.suffix.lower() == ".png":
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        for token in FORBIDDEN:
            if token in text:
                errors.append(f"private runtime token in {path.relative_to(atlas)}: {token}")
    return errors


def main() -> int:
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[2]
    errors = validate(root)
    print(json.dumps({"status": "PASS" if not errors else "FAIL", "errors": errors}, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
