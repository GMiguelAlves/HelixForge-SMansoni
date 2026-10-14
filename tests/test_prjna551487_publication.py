from __future__ import annotations

import hashlib
import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "results" / "rnaseq" / "PRJNA551487"


class PRJNA551487PublicationTest(unittest.TestCase):
    def test_public_core_and_scientific_exclusion(self) -> None:
        for relative in (
            "README.md",
            "qc/qc_sample_summary.tsv",
            "qc/qc_summary.json",
            "expression/counts_matrix.tsv",
            "expression/tpm_matrix.tsv",
            "expression/length_matrix.tsv",
            "expression/import_manifest.json",
            "manifests/rnaseq_run_manifest.json",
            "manifests/terminal_manifest_validation.json",
            "reports/PRJNA551487_multiqc.html",
        ):
            self.assertTrue((ROOT / relative).is_file(), relative)
        publication = json.loads((ROOT / "manifests/publication_validation.json").read_text(encoding="utf-8"))
        self.assertEqual("review_package", publication["status"])
        self.assertEqual("WITHHELD_AFTER_REVIEW", publication["de_promotion"])
        self.assertEqual("EXCLUDE_FROM_ATLAS", publication["atlas_inclusion"])
        self.assertEqual("10/10", publication["checks"]["qc_samples_review"])
        self.assertFalse((ROOT / "differential_expression").exists())

    def test_fullfile_diagnostic_denominators_and_warning(self) -> None:
        path = ROOT / "reports/fullfile_locus_qc_summary.tsv"
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(10, len(rows))
        self.assertEqual(10, len({row["sample_id"] for row in rows}))
        self.assertEqual(10, len({row["run_accession"] for row in rows}))
        self.assertEqual({"DMSO", "FPL"}, {row["arm"] for row in rows})
        self.assertEqual([5, 5], [sum(row["arm"] == arm for row in rows) for arm in ("DMSO", "FPL")])
        for row in rows:
            self.assertEqual(3000, int(row["random_pairs"]))
            self.assertGreater(int(row["full_fastq_pairs"]), 3000)
            self.assertGreater(int(row["sample_decile_min"]), 0)
            self.assertLessEqual(int(row["sample_decile_max"]), 3000)
            self.assertLess(float(row["salmon_mapping_percent_full_run"]), 60)
            self.assertLessEqual(int(row["locus_either_mate_hit"]), 3000)
            self.assertAlmostEqual(
                int(row["locus_either_mate_hit"]) / 3000,
                float(row["locus_either_mate_fraction"]),
            )
        report = (ROOT / "reports/PRJNA551487_execution_report.md").read_text(encoding="utf-8")
        self.assertIn("proof of a library-preparation mechanism", (ROOT / "README.md").read_text(encoding="utf-8"))
        self.assertIn("cannot be separated by this design", report)

    def test_import_matrices_match_native_manifest(self) -> None:
        manifest = json.loads((ROOT / "expression/import_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(10, manifest["sample_count"])
        for role in ("counts", "abundance", "lengths"):
            artifact = manifest["artifacts"][role]
            path = ROOT / "expression" / artifact["path"]
            self.assertEqual(artifact["sha256"], hashlib.sha256(path.read_bytes()).hexdigest(), role)

    def test_every_public_file_is_hashed_and_sanitized(self) -> None:
        checksums = (ROOT / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
        listed = set()
        for line in checksums:
            digest, relative = line.split("  ", 1)
            self.assertNotIn(relative, listed)
            path = ROOT / relative
            self.assertTrue(path.is_file(), relative)
            self.assertEqual(digest, hashlib.sha256(path.read_bytes()).hexdigest(), relative)
            listed.add(relative)
            if path.suffix.lower() in {".md", ".json", ".html", ".tsv", ".txt", ".yml"}:
                content = path.read_text(encoding="utf-8")
                for token in ("/scratch/", "/home/", "ra236875", "srv-slurm", "C:\\Users\\"):
                    self.assertNotIn(token, content, relative)
        actual = {
            str(path.relative_to(ROOT)).replace("\\", "/")
            for path in ROOT.rglob("*")
            if path.is_file() and path.name != "SHA256SUMS"
        }
        self.assertEqual(actual, listed)


if __name__ == "__main__":
    unittest.main()
