from __future__ import annotations

import csv
import hashlib
import json
import unittest
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    "PRJNA312093": (20, 20, "SINGLE", "full", "ENABLED"),
    "PRJNA494090": (24, 24, "PAIRED", "full", "ENABLED"),
    "PRJNA704646": (24, 8, "PAIRED", "import", "BLOCKED"),
    "PRJNA1150591": (15, 15, "PAIRED", "full", "ENABLED"),
    "PRJNA1177852": (8, 8, "PAIRED", "import", "BLOCKED"),
    "PRJNA551487": (10, 10, "PAIRED", "full", "ENABLED"),
    "PRJNA428470": (6, 6, "PAIRED", "full", "ENABLED"),
}


def read_rows(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


class P1Preparation(unittest.TestCase):
    def test_static_preparation_contract(self) -> None:
        pending = read_rows(ROOT / "metadata/rnaseq/p1_pending_manifest.tsv", "\t")
        self.assertEqual({row["study"] for row in pending}, set(EXPECTED))
        global_runs: set[str] = set()
        for study, (run_count, sample_count, layout, mode, de_status) in EXPECTED.items():
            with self.subTest(study=study):
                runs = read_rows(ROOT / "metadata" / "rnaseq" / study / "runs.tsv", "\t")
                samples = read_rows(ROOT / "metadata" / "rnaseq" / study / "samples.tsv", "\t")
                metadata = read_rows(ROOT / "metadata" / "rnaseq" / study / "metadata.csv", ",")
                self.assertEqual((len(runs), len(samples), len(metadata)),
                                 (run_count, sample_count, run_count))
                run_ids = {row["run_accession"] for row in runs}
                self.assertEqual(len(run_ids), run_count)
                self.assertFalse(global_runs & run_ids)
                global_runs.update(run_ids)
                self.assertEqual(run_ids, {row["run_accession"] for row in metadata})
                self.assertEqual({row["sample_id"] for row in runs},
                                 {row["sample_id"] for row in samples})
                self.assertEqual({row["library_layout"] for row in runs}, {layout})
                self.assertEqual({row["library_layout"] for row in metadata}, {layout.lower()})
                technical_runs = Counter(row["sample_id"] for row in runs)
                for sample in samples:
                    self.assertRegex(sample["sample_id"], r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
                    self.assertEqual(technical_runs[sample["sample_id"]],
                                     int(sample["technical_runs"]))
                for row in runs:
                    self.assertTrue(all(row[f"fastq_1_{field}"] for field in ("url", "md5", "bytes")))
                    self.assertRegex(row["fastq_1_md5"], r"^[0-9a-f]{32}$")
                    r2 = [row[f"fastq_2_{field}"] for field in ("url", "md5", "bytes")]
                    if layout == "SINGLE":
                        self.assertFalse(any(r2))
                    else:
                        self.assertTrue(all(r2))
                        self.assertRegex(row["fastq_2_md5"], r"^[0-9a-f]{32}$")
                provenance = json.loads((ROOT / "provenance" / study / "input_manifest.json").read_text(encoding="utf-8"))
                snapshot = ROOT / "provenance" / study / "source_snapshots/ena_read_run.tsv"
                self.assertEqual(provenance["source_snapshot"]["sha256"],
                                 hashlib.sha256(snapshot.read_bytes()).hexdigest())
                self.assertEqual((provenance["run_mode"], provenance["differential_expression"]),
                                 (mode, de_status))
                settings = (ROOT / "config" / study / "user_settings.sh").read_text(encoding="utf-8")
                self.assertIn(f'export HF_RNASEQ_RUN_MODE="{mode}"', settings)
                if de_status == "BLOCKED":
                    blocked = json.loads((ROOT / "config" / study / "de_spec.blocked.json").read_text(encoding="utf-8"))
                    self.assertEqual(blocked["status"], "BLOCKED")
                else:
                    spec = json.loads((ROOT / "config" / study / "de_spec.json").read_text(encoding="utf-8"))
                    levels = Counter(sample["condition"] for sample in samples)
                    for contrast in spec["contrasts"]:
                        for level in (contrast["numerator"], contrast["denominator"]):
                            self.assertGreaterEqual(levels[level], spec["parameters"]["min_replicates"])

    def test_single_end_prior_is_explicit(self) -> None:
        with (ROOT / "metadata/rnaseq/PRJNA312093/metadata.csv").open(encoding="utf-8", newline="") as handle:
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
