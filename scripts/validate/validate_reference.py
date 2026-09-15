#!/usr/bin/env python3
"""Validate the frozen S. mansoni reference without changing identifiers."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


VALID_DNA = set("ACGTURYSWKMBDHVNacgturyswkmbdhvn.-")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_stats(path: Path) -> tuple[dict[str, object], set[str]]:
    ids: list[str] = []
    total_bases = 0
    invalid = Counter()
    current = None
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                current = line[1:].split()[0]
                if not current:
                    raise ValueError(f"{path}: empty FASTA identifier at line {line_no}")
                ids.append(current)
            elif line:
                if current is None:
                    raise ValueError(f"{path}: sequence before first header")
                total_bases += len(line)
                invalid.update(char for char in line if char not in VALID_DNA)
    duplicates = sorted(key for key, count in Counter(ids).items() if count > 1)
    return {
        "path": path.name,
        "sha256": sha256(path),
        "sequence_count": len(ids),
        "total_bases": total_bases,
        "duplicate_ids": duplicates,
        "invalid_characters": dict(invalid),
    }, set(ids)


def parse_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for key, quoted, plain in re.findall(r'(\S+)\s+(?:"([^"]*)"|([^;\s]+))\s*;?', text):
        result[key] = quoted or plain
    return result


def annotation_stats(
    path: Path, gtf: bool
) -> tuple[dict[str, object], set[str], set[str], set[str]]:
    contigs: set[str] = set()
    genes: set[str] = set()
    transcripts: set[str] = set()
    features: Counter[str] = Counter()
    records = 0
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}: expected 9 columns at line {line_no}")
            records += 1
            contigs.add(fields[0])
            features[fields[2].lower()] += 1
            if gtf:
                attrs = parse_attributes(fields[8])
                if attrs.get("gene_id"):
                    genes.add(attrs["gene_id"])
                if attrs.get("transcript_id"):
                    transcripts.add(attrs["transcript_id"])
            else:
                attrs = dict(
                    item.split("=", 1) for item in fields[8].split(";") if "=" in item
                )
                kind = fields[2].lower()
                identifier = attrs.get("ID", "")
                if kind == "gene" and identifier:
                    genes.add(identifier.removeprefix("gene:"))
                if kind in {"mrna", "transcript"} and identifier:
                    transcripts.add(identifier.removeprefix("transcript:"))
    return {
        "path": path.name,
        "sha256": sha256(path),
        "records": records,
        "contig_count": len(contigs),
        "feature_counts": dict(sorted(features.items())),
        "gene_id_count": len(genes),
        "transcript_id_count": len(transcripts),
    }, contigs, genes, transcripts


def tx2gene_stats(path: Path, transcript_fasta_ids: set[str]) -> dict[str, object]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or list(rows[0]) != ["transcript_id", "gene_id"]:
        raise ValueError("tx2gene columns must be transcript_id and gene_id")
    pairs = [(row["transcript_id"], row["gene_id"]) for row in rows]
    if any(not transcript or not gene for transcript, gene in pairs):
        raise ValueError("tx2gene contains empty identifiers")
    by_transcript: dict[str, set[str]] = {}
    for transcript, gene in pairs:
        by_transcript.setdefault(transcript, set()).add(gene)
    conflicts = sorted(tx for tx, genes in by_transcript.items() if len(genes) > 1)
    mapped = set(by_transcript)
    return {
        "path": path.name,
        "sha256": sha256(path),
        "rows": len(rows),
        "genes": len({gene for _, gene in pairs}),
        "transcripts": len(mapped),
        "duplicate_pairs": len(pairs) - len(set(pairs)),
        "multi_gene_transcripts": conflicts,
        "transcriptome_mapped": len(transcript_fasta_ids & mapped),
        "transcriptome_unmapped": len(transcript_fasta_ids - mapped),
        "annotation_only_transcripts": len(mapped - transcript_fasta_ids),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = args.reference_dir.resolve()
    genome = root / "schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa"
    transcriptome = root / "schistosoma_mansoni.PRJEA36577.WBPS19.mRNA_transcripts.fa"
    gff3 = root / "schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3"
    gtf = root / "schistosoma_mansoni.PRJEA36577.WBPS19.canonical_geneset.gtf"
    tx2gene = root / "tx2gene.tsv"

    genome_result, genome_ids = fasta_stats(genome)
    transcript_result, transcript_ids = fasta_stats(transcriptome)
    gff_result, gff_contigs, _, _ = annotation_stats(gff3, gtf=False)
    gtf_result, gtf_contigs, _, _ = annotation_stats(gtf, gtf=True)
    mapping_result = tx2gene_stats(tx2gene, transcript_ids)
    missing_gff = sorted(gff_contigs - genome_ids)
    missing_gtf = sorted(gtf_contigs - genome_ids)
    index_info = json.loads((root / "salmon_index_k31/versionInfo.json").read_text())

    errors = []
    for label, result in (("genome", genome_result), ("transcriptome", transcript_result)):
        if not result["sequence_count"] or result["duplicate_ids"] or result["invalid_characters"]:
            errors.append(f"{label} FASTA sanity failed")
    if not gff_result["records"] or not gtf_result["records"]:
        errors.append("annotation is empty")
    if missing_gff or missing_gtf:
        errors.append("annotation contains contigs absent from genome FASTA")
    if mapping_result["duplicate_pairs"] or mapping_result["multi_gene_transcripts"]:
        errors.append("tx2gene mapping is not one-to-one by transcript")
    if mapping_result["transcriptome_unmapped"]:
        errors.append("official transcriptome contains unmapped transcripts")
    if index_info.get("salmonVersion") != "1.10.3" or index_info.get("auxKmerLength") != 31:
        errors.append("Salmon index version or k-mer size mismatch")

    result = {
        "schema_version": "1.0",
        "status": "PASS" if not errors else "FAIL",
        "reference_id": "Schistosoma_mansoni_SM_V10_WBPS19",
        "genome": genome_result,
        "transcriptome": transcript_result,
        "gff3": gff_result,
        "gtf": gtf_result,
        "contig_compatibility": {
            "status": "PASS" if not missing_gff and not missing_gtf else "FAIL",
            "gff3_missing_from_genome": missing_gff,
            "gtf_missing_from_genome": missing_gtf,
        },
        "tx2gene": mapping_result,
        "salmon_index": index_info,
        "errors": errors,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
