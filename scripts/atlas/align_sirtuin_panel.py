"""Align a focused sirtuin panel and trim sparsely occupied columns."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


def fasta(path: Path) -> dict[str, str]:
    result = {}
    name, chunks = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name:
                result[name] = "".join(chunks)
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if name:
        result[name] = "".join(chunks)
    return result


def write_fasta(path: Path, data: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in data.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--vendor-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.vendor_dir.resolve()))
    import pyfamsa  # noqa: PLC0415

    args.output.mkdir(parents=True, exist_ok=True)
    source = fasta(args.input)
    msa = pyfamsa.Aligner().align([pyfamsa.Sequence(name.encode(), sequence.encode()) for name, sequence in source.items()])
    aligned = {record.id.decode(): record.sequence.decode() for record in msa}
    lengths = {len(sequence) for sequence in aligned.values()}
    if len(lengths) != 1:
        raise ValueError(f"unequal aligned lengths: {lengths}")
    untrimmed = next(iter(lengths))
    write_fasta(args.output / "sirtuins_aligned.afa", aligned)
    retain = [index for index in range(untrimmed)
              if sum(sequence[index] != "-" for sequence in aligned.values()) / len(aligned) >= 0.7]
    if len(retain) < 150:
        raise ValueError(f"only {len(retain)} columns at >=70% occupancy")
    trimmed = {name: "".join(sequence[i] for i in retain) for name, sequence in aligned.items()}
    write_fasta(args.output / "sirtuins_70pct_occupancy.afa", trimmed)
    summary = {"method": "pyfamsa", "version": pyfamsa.__version__, "input_sequences": len(source),
               "aligned_columns": untrimmed, "retained_columns": len(retain), "occupancy_threshold": 0.7,
               "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
               "trimmed_sha256": hashlib.sha256((args.output / "sirtuins_70pct_occupancy.afa").read_bytes()).hexdigest()}
    (args.output / "alignment_manifest.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
