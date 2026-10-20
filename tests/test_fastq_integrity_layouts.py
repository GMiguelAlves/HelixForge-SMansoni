from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIELDS = (
    "dataset", "sample_id", "run_accession", "library_layout",
    "fastq_1_url", "fastq_2_url", "fastq_1_md5", "fastq_2_md5",
    "fastq_1_bytes", "fastq_2_bytes",
)


class FastqIntegrityLayouts(unittest.TestCase):
    def check_layout(self, layout: str) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            package = root / "package"
            data = root / "data"
            manifest_dir = package / "metadata" / "rnaseq" / "TEST"
            fastq_dir = data / "TEST" / "fastq_ftp"
            manifest_dir.mkdir(parents=True)
            fastq_dir.mkdir(parents=True)
            row = dict.fromkeys(FIELDS, "")
            row.update(dataset="TEST", sample_id="sample1", run_accession="RUN1", library_layout=layout)
            for mate in (("1",) if layout == "SINGLE" else ("1", "2")):
                content = gzip.compress(f"@read{mate}\nAC\n+\n!!\n".encode())
                (fastq_dir / f"RUN1_{mate}.fastq.gz").write_bytes(content)
                row[f"fastq_{mate}_url"] = f"https://example.org/RUN1_{mate}.fastq.gz"
                row[f"fastq_{mate}_md5"] = hashlib.md5(content).hexdigest()
                row[f"fastq_{mate}_bytes"] = str(len(content))
            with (manifest_dir / "runs.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
                writer.writeheader()
                writer.writerow(row)

            output = root / "audit"
            result = subprocess.run(
                [sys.executable, str(ROOT / "scripts/validate/validate_fastq_integrity.py"),
                 "TEST", "--package-root", str(package), "--data-root", str(data),
                 "--output-json", str(output / "integrity.json"),
                 "--output-tsv", str(output / "integrity.tsv"),
                 "--output-sha256", str(output / "SHA256SUMS")],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            report = json.loads((output / "integrity.json").read_text(encoding="utf-8"))
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(report["layout"], layout)
            self.assertEqual(report["fastq_files_expected"], 1 if layout == "SINGLE" else 2)
            self.assertEqual(report["r2_total_bytes"], 0 if layout == "SINGLE" else int(row["fastq_2_bytes"]))

            legacy = subprocess.run(
                [sys.executable, str(ROOT / "scripts/validate/validate_downloads.py"),
                 "TEST", "--package-root", str(package), "--data-root", str(data)],
                capture_output=True, text=True,
            )
            self.assertEqual(legacy.returncode, 0, legacy.stdout + legacy.stderr)

    def test_single_end(self) -> None:
        self.check_layout("SINGLE")

    def test_paired_end(self) -> None:
        self.check_layout("PAIRED")


if __name__ == "__main__":
    unittest.main()
