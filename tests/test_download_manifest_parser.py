from __future__ import annotations

import csv
import shutil
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "scripts/download/parse_manifest_row.sh"
PROJECTS = (
    "PRJNA312093",
    "PRJNA494090",
    "PRJNA704646",
    "PRJNA1150591",
    "PRJNA1177852",
    "PRJNA551487",
    "PRJNA428470",
)
FIELDS = (
    "dataset", "sample_id", "run_accession", "biosample", "secondary_biosample",
    "experiment_accession", "technical_unit", "library_name", "library_layout",
    "fastq_1_url", "fastq_2_url", "fastq_1_md5", "fastq_2_md5",
    "fastq_1_bytes", "fastq_2_bytes",
)


@unittest.skipUnless(sys.platform != "win32" and shutil.which("bash"), "requires Bash")
class DownloadManifestParser(unittest.TestCase):
    def test_first_and_last_run_of_each_p1_study(self) -> None:
        for project in PROJECTS:
            manifest = ROOT / "metadata" / "rnaseq" / project / "runs.tsv"
            with manifest.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(rows, project)
            for task in {1, len(rows)}:
                with self.subTest(project=project, task=task):
                    command = (
                        'source "$1"; download_parse_manifest_row "$2" "$3"; '
                        'printf "%s\\0" "$dataset" "$sample" "$run" "$biosample" '
                        '"$secondary" "$experiment" "$unit" "$library" "$layout" '
                        '"$url1" "$url2" "$md51" "$md52" "$bytes1" "$bytes2"'
                    )
                    result = subprocess.run(
                        ["bash", "-c", command, "bash", str(HELPER), str(manifest), str(task)],
                        check=True,
                        capture_output=True,
                    )
                    actual = result.stdout.decode("utf-8").split("\0")[:-1]
                    expected = [rows[task - 1][field] for field in FIELDS]
                    self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
