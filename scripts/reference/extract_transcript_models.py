#!/usr/bin/env python3
"""Extract exon structures from the checksum-verified frozen WBPS GFF3."""

import argparse
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote


def attributes(text):
    return {key: unquote(value) for item in text.rstrip(";").split(";")
            if "=" in item for key, value in [item.split("=", 1)]}


def extract_models(path, expected_sha256):
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected_sha256:
        raise ValueError("GFF3 checksum differs from the frozen reference manifest")
    transcripts = {}
    exons = defaultdict(list)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"mRNA", "exon"}:
                continue
            attr = attributes(fields[8])
            start, end = int(fields[3]), int(fields[4])
            if fields[2] == "mRNA":
                transcript = attr["ID"]
                if transcript in transcripts:
                    raise ValueError(f"duplicate transcript: {transcript}")
                transcripts[transcript] = {
                    "gene_id": attr["Parent"].removeprefix("gene:"),
                    "transcript_id": transcript.removeprefix("transcript:"),
                    "chromosome": fields[0], "start": start, "end": end,
                    "strand": fields[6],
                }
            else:
                for parent in attr["Parent"].split(","):
                    exons[parent].append((fields[0], start, end, fields[6]))
    genes = defaultdict(list)
    for transcript, model in sorted(transcripts.items()):
        intervals = sorted(set(exons[transcript]), key=lambda row: (row[1], row[2]))
        if not intervals:
            raise ValueError(f"no annotated exons: {transcript}")
        previous_end = 0
        for chrom, start, end, strand in intervals:
            if chrom != model["chromosome"] or strand != model["strand"]:
                raise ValueError(f"exon identity differs from transcript: {transcript}")
            if not model["start"] <= start <= end <= model["end"] or start <= previous_end:
                raise ValueError(f"invalid or overlapping exons: {transcript}")
            previous_end = end
        gene = model.pop("gene_id")
        model["exons"] = [[row[1], row[2]] for row in intervals]
        genes[gene].append(model)
    return dict(sorted(genes.items()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    reference = json.loads(args.reference_manifest.read_text(encoding="utf-8"))
    source = next(row for row in reference["artifacts"] if row["role"] == "annotation_gff3")
    genes = extract_models(args.gff, source["sha256"])
    result = {
        "schema_version": "1.0", "reference_id": reference["reference_id"],
        "coordinate_system": "1-based-inclusive",
        "source_url": source["url"], "source_sha256": source["sha256"],
        "display_policy": "all_annotated_transcripts_separate_rows",
        "genes": genes,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(",", ":")) + "\n", encoding="utf-8")
    print(json.dumps({"genes": len(genes), "transcripts": sum(map(len, genes.values())),
                      "source_sha256": source["sha256"]}))


if __name__ == "__main__":
    main()
