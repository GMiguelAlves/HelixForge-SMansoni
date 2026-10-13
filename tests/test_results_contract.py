from __future__ import annotations

import unittest
import json
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "rnaseq"
STUDIES = (
    "PRJNA602528",
    "PRJNA597909",
    "PRJEB14695",
    "PRJEB32839",
    "E-MTAB-451",
    "PRJEB3190",
    "E-ERAD-478",
    "PRJNA312093",
    "PRJNA1150591",
    "PRJNA1177852",
    "PRJNA704646",
)


class PublicResultsContract(unittest.TestCase):
    def test_all_studies_expose_the_stable_core(self) -> None:
        required = (
            "README.md",
            "qc/qc_sample_summary.tsv",
            "qc/qc_summary.json",
            "expression/counts_matrix.tsv",
            "expression/tpm_matrix.tsv",
            "expression/length_matrix.tsv",
        )
        for study in STUDIES:
            with self.subTest(study=study):
                for relative in required:
                    self.assertTrue(
                        (RESULTS / study / relative).is_file(),
                        f"{study}: missing stable result {relative}",
                    )
                if study in {"PRJNA1177852", "PRJNA704646"}:
                    self.assertTrue((RESULTS / study / "expression/import_manifest.json").is_file())
                    self.assertTrue((RESULTS / study / "manifests/publication_validation.json").is_file())
                    self.assertFalse((RESULTS / study / "manifests/rnaseq_run_manifest.json").exists())
                else:
                    self.assertTrue((RESULTS / study / "manifests/rnaseq_run_manifest.json").is_file())
                    self.assertTrue((RESULTS / study / "manifests/terminal_manifest_validation.json").is_file())

    def test_e_erad_478_uses_the_compact_public_layout(self) -> None:
        root = RESULTS / "E-ERAD-478"
        self.assertFalse((root / "compatibility").exists())
        self.assertFalse((root / "rnaseq").exists())
        self.assertTrue((root / "reports/E-ERAD-478_multiqc.html").is_file())
        self.assertTrue((root / "reports/E-ERAD-478_execution_report.md").is_file())
        self.assertTrue((root / "reports/final_validation.json").is_file())

    def test_multiqc_is_not_mixed_with_the_compact_qc_interface(self) -> None:
        root = RESULTS / "E-ERAD-478"
        self.assertEqual([], list((root / "qc").glob("*_multiqc.html")))

    def test_public_e_erad_has_no_large_reentry_tree(self) -> None:
        root = RESULTS / "E-ERAD-478"
        forbidden_parts = {"integration_artifacts", "quants"}
        offenders = [
            str(path.relative_to(root))
            for path in root.rglob("*")
            if path.is_file() and forbidden_parts.intersection(path.parts)
        ]
        self.assertEqual([], offenders)

    def test_public_studies_exclude_salmon_reentry_payloads(self) -> None:
        offenders = []
        for study in STUDIES:
            root = RESULTS / study
            offenders.extend(
                str(path.relative_to(RESULTS))
                for path in root.rglob("*")
                if path.is_file()
                and (
                    path.name == "quant.sf"
                    or "quants" in path.relative_to(root).parts
                )
            )
        self.assertEqual([], offenders)

    def test_prjna312093_remains_on_biological_hold(self) -> None:
        root = RESULTS / "PRJNA312093"
        gate = json.loads((root / "reports/final_validation.json").read_text(encoding="utf-8"))
        self.assertEqual("PASS", gate["technical_qc"])
        self.assertEqual("REVIEW", gate["biological_replicate_concordance"])
        self.assertEqual("HOLD", gate["study_promotion"])
        self.assertEqual("HOLD", gate["atlas_inclusion"])
        self.assertTrue((root / "reports/PRJNA312093_multiqc.html").is_file())
        self.assertFalse((root / "differential_expression/dds_condition.rds").exists())

    def test_prjna312093_public_files_are_sanitized_and_hashed(self) -> None:
        root = RESULTS / "PRJNA312093"
        entries = (root / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
        listed = set()
        for line in entries:
            digest, relative = line.split("  ", 1)
            path = root / relative
            self.assertTrue(path.is_file(), relative)
            self.assertEqual(digest, hashlib.sha256(path.read_bytes()).hexdigest())
            listed.add(relative)
        actual = {str(path.relative_to(root)).replace("\\", "/") for path in root.rglob("*") if path.is_file() and path.name != "SHA256SUMS"}
        self.assertEqual(actual, listed)
        forbidden = ("/scratch/", "/home/", "ra236875", "srv-slurm", "C:\\Users\\")
        for relative in sorted(actual):
            path = root / relative
            if path.suffix.lower() in {".md", ".json", ".html", ".tsv", ".txt", ".yml"}:
                text = path.read_text(encoding="utf-8")
                for token in forbidden:
                    self.assertNotIn(token, text, relative)

    def test_prjna312093_matrices_match_native_import_manifest(self) -> None:
        root = RESULTS / "PRJNA312093"
        manifest = json.loads((root / "expression/import_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(20, manifest["sample_count"])
        for role in ("counts", "abundance", "lengths"):
            artifact = manifest["artifacts"][role]
            path = root / "expression" / artifact["path"]
            self.assertEqual(artifact["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
        de = json.loads((root / "differential_expression/de_manifest.json").read_text(encoding="utf-8"))
        significant = de["artifacts"]["significant"]
        path = root / "differential_expression" / significant["path"]
        self.assertEqual(significant["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())

    def test_prjna1150591_public_package(self) -> None:
        root = RESULTS / "PRJNA1150591"
        hashes = dict(line.split("  ", 1) for line in (root / "SHA256SUMS").read_text(encoding="utf-8").splitlines())
        actual = {
            str(path.relative_to(root)).replace("\\", "/")
            for path in root.rglob("*")
            if path.is_file() and path.name != "SHA256SUMS"
        }
        self.assertEqual(actual, set(hashes.values()))
        for digest, relative in hashes.items():
            path = root / relative
            self.assertEqual(digest, hashlib.sha256(path.read_bytes()).hexdigest(), relative)
            if path.suffix.lower() in {".md", ".json", ".html", ".tsv", ".txt", ".yml"}:
                content = path.read_text(encoding="utf-8")
                for forbidden in ("/scratch/", "/home/", "ra236875", "srv-slurm", "C:\\Users\\"):
                    self.assertNotIn(forbidden, content, relative)
        import_manifest = json.loads((root / "expression/import_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(15, import_manifest["sample_count"])
        for role in ("counts", "abundance", "lengths"):
            artifact = import_manifest["artifacts"][role]
            path = root / "expression" / artifact["path"]
            self.assertEqual(artifact["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
        de = json.loads((root / "differential_expression/de_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual("~ condition", de["design"]["formula"])
        self.assertEqual(1, de["contrasts"])
        self.assertEqual("complete", de["status"])

    def test_prjna1177852_import_only_package(self) -> None:
        root = RESULTS / "PRJNA1177852"
        hashes = dict(line.split("  ", 1) for line in (root / "SHA256SUMS").read_text(encoding="utf-8").splitlines())
        actual = {
            str(path.relative_to(root)).replace("\\", "/")
            for path in root.rglob("*")
            if path.is_file() and path.name != "SHA256SUMS"
        }
        self.assertEqual(actual, set(hashes.values()))
        for digest, relative in hashes.items():
            path = root / relative
            self.assertEqual(digest, hashlib.sha256(path.read_bytes()).hexdigest(), relative)
            if path.suffix.lower() in {".md", ".json", ".html", ".tsv", ".txt", ".yml"}:
                content = path.read_text(encoding="utf-8")
                for forbidden in ("/scratch/", "/home/", "ra236875", "srv-slurm", "C:\\Users\\"):
                    self.assertNotIn(forbidden, content, relative)
        manifest = json.loads((root / "expression/import_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(8, manifest["sample_count"])
        for role in ("counts", "abundance", "lengths"):
            artifact = manifest["artifacts"][role]
            self.assertEqual(
                artifact["sha256"],
                hashlib.sha256((root / "expression" / artifact["path"]).read_bytes()).hexdigest(),
            )
        publication = json.loads((root / "manifests/publication_validation.json").read_text(encoding="utf-8"))
        self.assertEqual("import", publication["run_mode"])
        self.assertEqual("DESCRIPTIVE_ONLY", publication["scientific_classification"])
        self.assertEqual("not_run_by_design", publication["differential_expression"])
        self.assertFalse((root / "differential_expression").exists())

    def test_prjna704646_descriptive_package_keeps_qc_review(self) -> None:
        root = RESULTS / "PRJNA704646"
        hashes = dict(line.split("  ", 1) for line in (root / "SHA256SUMS").read_text(encoding="utf-8").splitlines())
        actual = {
            str(path.relative_to(root)).replace("\\", "/")
            for path in root.rglob("*")
            if path.is_file() and path.name != "SHA256SUMS"
        }
        self.assertEqual(actual, set(hashes.values()))
        for digest, relative in hashes.items():
            path = root / relative
            self.assertEqual(digest, hashlib.sha256(path.read_bytes()).hexdigest(), relative)
            if path.suffix.lower() in {".md", ".json", ".html", ".tsv", ".txt", ".yml"}:
                content = path.read_text(encoding="utf-8")
                for forbidden in ("/scratch/", "/home/", "ra236875", "srv-slurm", "C:\\Users\\"):
                    self.assertNotIn(forbidden, content, relative)
        manifest = json.loads((root / "expression/import_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(8, manifest["sample_count"])
        for role in ("counts", "abundance", "lengths"):
            artifact = manifest["artifacts"][role]
            self.assertEqual(
                artifact["sha256"],
                hashlib.sha256((root / "expression" / artifact["path"]).read_bytes()).hexdigest(),
            )
        qc = json.loads((root / "qc/qc_summary.json").read_text(encoding="utf-8"))
        self.assertEqual({"PASS": 5, "REVIEW": 3}, qc["sample_classification"])
        publication = json.loads((root / "manifests/publication_validation.json").read_text(encoding="utf-8"))
        self.assertEqual("DESCRIPTIVE_ONLY", publication["scientific_classification"])
        self.assertEqual("not_run_by_design", publication["differential_expression"])
        self.assertFalse((root / "differential_expression").exists())


if __name__ == "__main__":
    unittest.main()
