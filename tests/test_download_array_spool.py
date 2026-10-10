from __future__ import annotations

import csv
import gzip
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIELDS = (
    "dataset", "sample_id", "run_accession", "biosample", "secondary_biosample",
    "experiment_accession", "technical_unit", "library_name", "library_layout",
    "fastq_1_url", "fastq_2_url", "fastq_1_md5", "fastq_2_md5",
    "fastq_1_bytes", "fastq_2_bytes",
)


@unittest.skipUnless(sys.platform != "win32" and shutil.which("bash"), "requires Bash")
class DownloadArraySpool(unittest.TestCase):
    def test_single_end_script_runs_from_slurm_spool_location(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            package = root / "package"
            source_dir = package / "scripts/download"
            source_dir.mkdir(parents=True)
            shutil.copyfile(ROOT / "scripts/download/parse_manifest_row.sh", source_dir / "parse_manifest_row.sh")
            spool = root / "spool"
            spool.mkdir()
            shutil.copyfile(ROOT / "scripts/download/download_fastq_array.sh", spool / "slurm_script")
            fixture = root / "source.fastq.gz"
            fixture.write_bytes(gzip.compress(b"@read\nAC\n+\n!!\n"))
            row = dict.fromkeys(FIELDS, "")
            row.update(
                dataset="TEST", sample_id="sample1", run_accession="RUN1", library_layout="SINGLE",
                fastq_1_url=str(fixture), fastq_1_md5=hashlib.md5(fixture.read_bytes()).hexdigest(),
                fastq_1_bytes=str(fixture.stat().st_size),
            )
            manifest = root / "runs.tsv"
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
                writer.writeheader()
                writer.writerow(row)
            fake_bin = root / "bin"
            fake_bin.mkdir()
            curl = fake_bin / "curl"
            curl.write_text(
                '#!/usr/bin/env bash\nwhile [[ $# -gt 0 ]]; do '
                'if [[ "$1" == -o ]]; then output=$2; shift 2; '
                'else source_file=$1; shift; fi; done\ncp "$source_file" "$output"\n',
                encoding="utf-8",
            )
            curl.chmod(0o755)
            data = root / "data"
            env = os.environ | {
                "HF_PACKAGE_ROOT": str(package), "HF_DATA_ROOT": str(data),
                "SLURM_ARRAY_TASK_ID": "1", "PATH": f"{fake_bin}:{os.environ['PATH']}",
            }
            result = subprocess.run(["bash", str(spool / "slurm_script"), str(manifest)],
                                    env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual((data / "TEST/fastq_ftp/RUN1_1.fastq.gz").read_bytes(), fixture.read_bytes())
            self.assertTrue((data / "TEST/download_provenance/RUN1.json").is_file())


if __name__ == "__main__":
    unittest.main()
