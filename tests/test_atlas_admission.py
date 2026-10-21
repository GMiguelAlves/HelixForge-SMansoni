import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts/atlas"))
from build_rnaseq_atlas import approved_contrast_sources, select_expression_samples, study_spec


class AdmissionTests(unittest.TestCase):
    def test_frozen_sample_selection_preserves_values_and_order(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as directory:
            path = Path(directory) / "samples.tsv"
            path.write_text("sample_key\nTEST__female\nTEST__male\n", encoding="utf-8")
            columns, values = select_expression_samples(
                ["TEST__adult", "TEST__male", "TEST__female"],
                np.array([[99, 2, 3], [88, 4, 5]]), path, "TEST"
            )
            self.assertEqual(columns, ["TEST__female", "TEST__male"])
            np.testing.assert_array_equal(values, [[3, 2], [5, 4]])

    def test_missing_or_duplicate_approved_sample_is_rejected(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as directory:
            path = Path(directory) / "samples.tsv"
            for content in ("sample_key\nTEST__missing\n", "sample_key\nTEST__a\nTEST__a\n"):
                path.write_text(content, encoding="utf-8")
                with self.assertRaises(ValueError):
                    select_expression_samples(["TEST__a"], np.array([[1]]), path, "TEST")

    def test_approved_refit_hash_failure_does_not_fall_back(self):
        entry = json.loads((ROOT / "analyses/atlas/atlas_config.json").read_text(encoding="utf-8"))["studies"]
        entry = next(row for row in entry if row["id"] == "PRJNA312093")
        _, spec = study_spec(ROOT, entry)
        self.assertEqual(len(spec["contrasts"]), 3)
        sources = approved_contrast_sources(ROOT, entry, spec)
        self.assertTrue(all("restricted_refit" in str(path) for _, path in sources))
        entry["contrast_sources"][spec["contrasts"][0]["id"]]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            approved_contrast_sources(ROOT, entry, spec)

    def test_unknown_approved_contrast_is_rejected(self):
        entry = {"id": "PRJNA312093", "differential_expression": True,
                 "contrast_sources": {"not_a_frozen_contrast": {}}, "analysis_id": "test"}
        with self.assertRaises(ValueError):
            study_spec(ROOT, entry)


if __name__ == "__main__":
    unittest.main()
