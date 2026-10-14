from __future__ import annotations

import csv
import hashlib
import json
import unittest
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1] / "results" / "rnaseq" / "PRJNA312093"
REFIT = PACKAGE / "restricted_refit"


class RestrictedRefitPublicationTest(unittest.TestCase):
    def test_only_reviewed_stage_contrasts_are_approved(self) -> None:
        manifest = json.loads((REFIT / "manifest.json").read_text(encoding="utf-8"))
        self.assertEqual("PRJNA312093", manifest["study"])
        self.assertEqual("atlas_approved_restricted_refit", manifest["type"])
        self.assertEqual("~ condition", manifest["design"])
        self.assertEqual(12, manifest["sample_count"])
        self.assertEqual("included_samples.tsv", manifest["approved_expression_samples"])
        self.assertEqual("select_only_the_12_included_samples", manifest["expression_rule"])
        self.assertEqual(9342, manifest["genes_tested_per_contrast"])
        self.assertEqual([1, 2, 3], [item["stage"] for item in manifest["contrasts"]])
        self.assertEqual([290, 1045, 437], [item["significant_genes"] for item in manifest["contrasts"]])
        self.assertEqual(
            {"RESTRICTED_REFIT", "TWO_REPLICATES_PER_GROUP", "SOURCE_IDENTITY_LIMITATION"},
            set(manifest["warnings"]),
        )
        self.assertEqual(3, len(list((REFIT / "contrasts").glob("*.tsv"))))
        for item in manifest["contrasts"]:
            path = REFIT / item["path"]
            self.assertEqual(item["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
            with path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(9342, len(rows))
            significant = sum(
                bool(row["padj"])
                and float(row["padj"]) < 0.05
                and abs(float(row["log2FoldChange"])) >= 1
                for row in rows
            )
            self.assertEqual(item["significant_genes"], significant)

    def test_sample_contract_and_audit_boundary(self) -> None:
        with (REFIT / "included_samples.tsv").open(encoding="utf-8", newline="") as handle:
            samples = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(12, len(samples))
        self.assertEqual(
            {(f"schistosomulum_stage_{stage}", sex) for stage in range(1, 4) for sex in ("female", "male")},
            {(row["stage"], row["sex"]) for row in samples},
        )
        self.assertTrue(all(row["stage"].startswith("schistosomulum_stage_") for row in samples))
        with (PACKAGE / "expression/tpm_matrix.tsv").open(encoding="utf-8") as handle:
            header = handle.readline().rstrip("\n\r").split("\t")
        self.assertEqual("gene_id", header[0])
        self.assertTrue({row["__rowname"] for row in samples} <= set(header[1:]))
        for path in REFIT.rglob("*"):
            if path.is_file() and path.suffix in {".md", ".json", ".tsv"}:
                value = path.read_text(encoding="utf-8")
                self.assertNotIn("/home/", value)
                self.assertNotIn("/scratch/", value)
                self.assertNotIn("C:\\Users\\", value)
        self.assertNotIn("adult", " ".join(item["path"] for item in json.loads(
            (REFIT / "manifest.json").read_text(encoding="utf-8")
        )["contrasts"]))

    def test_root_package_checksums_cover_refit(self) -> None:
        listed = {}
        for entry in (PACKAGE / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            digest, relative = entry.split("  ", 1)
            self.assertNotIn(relative, listed)
            self.assertEqual(digest, hashlib.sha256((PACKAGE / relative).read_bytes()).hexdigest(), relative)
            listed[relative] = digest
        actual = {
            path.relative_to(PACKAGE).as_posix()
            for path in PACKAGE.rglob("*")
            if path.is_file() and path.name != "SHA256SUMS"
        }
        self.assertEqual(actual, set(listed))


if __name__ == "__main__":
    unittest.main()
