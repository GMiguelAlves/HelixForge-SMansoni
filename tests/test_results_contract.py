from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
STUDIES = (
    "PRJNA602528",
    "PRJNA597909",
    "PRJEB14695",
    "PRJEB32839",
    "E-MTAB-451",
    "PRJEB3190",
    "E-ERAD-478",
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
            "manifests/rnaseq_run_manifest.json",
            "manifests/terminal_manifest_validation.json",
        )
        for study in STUDIES:
            with self.subTest(study=study):
                for relative in required:
                    self.assertTrue(
                        (RESULTS / study / relative).is_file(),
                        f"{study}: missing stable result {relative}",
                    )

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


if __name__ == "__main__":
    unittest.main()
