from __future__ import annotations

import csv
import hashlib
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "results" / "rnaseq" / "PRJNA428470"


class PRJNA428470PublicationTest(unittest.TestCase):
    def test_scientific_scope_and_sample_gate(self) -> None:
        qc = json.loads((ROOT / "qc/qc_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(6, qc["samples"])
        self.assertEqual({"PASS": 2, "REVIEW": 4}, qc["sample_classification"])
        self.assertEqual(60.0, qc["thresholds"]["salmon_mapping_review_below_percent"])
        self.assertLess(qc["salmon_mapping_percent_min"], 60.0)

        with (ROOT / "qc/qc_sample_summary.tsv").open(encoding="utf-8", newline="") as handle:
            samples = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(6, len(samples))
        self.assertEqual(4, sum(row["classification"] == "REVIEW" for row in samples))
        self.assertEqual(3, sum("5azac" in row["sample_id"] for row in samples))
        self.assertEqual(3, sum("control" in row["sample_id"] for row in samples))

        de = json.loads((ROOT / "differential_expression/de_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual("~ condition", de["design"]["formula"])
        self.assertEqual(1, de["contrasts"])
        self.assertEqual("complete", de["status"])
        with (ROOT / "differential_expression/contrasts/DEG_condition__female_adult_5azac_vs_female_adult_control.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(9307, len(rows))
        significant = sum(
            bool(row["padj"])
            and row["padj"] != "NA"
            and float(row["padj"]) < 0.05
            and abs(float(row["log2FoldChange"])) >= 1
            for row in rows
        )
        self.assertEqual(1454, significant)

    def test_supplementary_sensitivity_is_labelled_and_complete(self) -> None:
        with (ROOT / "qc/leave_one_out_summary.tsv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(6, len(rows))
        self.assertEqual(6, len({row["omitted_sample"] for row in rows}))
        for row in rows:
            self.assertEqual(5, int(row["retained_samples"]))
            self.assertEqual(1454, int(row["original_significant"]))
            self.assertGreaterEqual(int(row["original_significant_retained"]), 1286)
            self.assertEqual(1.0, float(row["original_significant_sign_agreement"]))
            self.assertGreaterEqual(float(row["all_gene_lfc_spearman"]), 0.976)
        review = (ROOT / "reports/PRJNA428470_execution_report.md").read_text(encoding="utf-8")
        self.assertIn("not a new primary model", review)
        self.assertIn("mRNA", review)
        self.assertIn("REVIEW_4_OF_6", review)

    def test_public_matrices_match_import_manifest(self) -> None:
        manifest = json.loads((ROOT / "expression/import_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(6, manifest["sample_count"])
        for role in ("counts", "abundance", "lengths"):
            artifact = manifest["artifacts"][role]
            path = ROOT / "expression" / artifact["path"]
            self.assertEqual(artifact["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())

    def test_compact_package_is_hashed_and_private_paths_are_absent(self) -> None:
        entries = (ROOT / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
        listed: set[str] = set()
        for entry in entries:
            digest, relative = entry.split("  ", 1)
            self.assertNotIn(relative, listed)
            path = ROOT / relative
            self.assertTrue(path.is_file(), relative)
            self.assertEqual(digest, hashlib.sha256(path.read_bytes()).hexdigest(), relative)
            listed.add(relative)
            if path.suffix.lower() in {".md", ".json", ".html", ".tsv", ".txt", ".yml"}:
                text = path.read_text(encoding="utf-8")
                for token in ("/scratch/", "/home/", "srv-slurm", "C:\\Users\\"):
                    self.assertNotIn(token, text, relative)
        actual = {
            str(path.relative_to(ROOT)).replace("\\", "/")
            for path in ROOT.rglob("*")
            if path.is_file() and path.name != "SHA256SUMS"
        }
        self.assertEqual(actual, listed)
        self.assertFalse((ROOT / "compatibility").exists())
        self.assertFalse((ROOT / "rnaseq").exists())
        self.assertFalse((ROOT / "differential_expression/dds_condition.rds").exists())
        self.assertTrue((ROOT / "reports/PRJNA428470_multiqc.html").is_file())
        self.assertEqual(
            "valid",
            json.loads((ROOT / "manifests/terminal_manifest_validation.json").read_text(encoding="utf-8"))["schema"],
        )


if __name__ == "__main__":
    unittest.main()
