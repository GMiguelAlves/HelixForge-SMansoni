from __future__ import annotations

import csv
import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "validate"))
from validate_p1_preparation import validate  # noqa: E402


class P1Preparation(unittest.TestCase):
    def test_static_preparation_contract(self) -> None:
        self.assertEqual(validate(ROOT), [])

    def test_single_end_prior_is_explicit(self) -> None:
        with (ROOT / "metadata/PRJNA312093/metadata.csv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 20)
        self.assertEqual({row["library_layout"] for row in rows}, {"single"})
        self.assertEqual({row["fragment_length_mean"] for row in rows}, {"200"})
        self.assertEqual({row["fragment_length_sd"] for row in rows}, {"80"})

    def test_frozen_gene_set_has_no_unresolved_historical_ids(self) -> None:
        manifest = json.loads((ROOT / "provenance/histone_modifier_gene_set_v2.json").read_text(encoding="utf-8"))
        mappings = manifest["counts"]["historical_mapping_associations"]
        self.assertEqual(manifest["status"], "FROZEN")
        self.assertNotIn("AMBIGUOUS_PREVIOUS_STABLE_ID", mappings)
        self.assertNotIn("NOT_FOUND_IN_WBPS19", mappings)


if __name__ == "__main__":
    unittest.main()
