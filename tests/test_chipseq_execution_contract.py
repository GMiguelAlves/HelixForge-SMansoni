from __future__ import annotations

import csv
import json
import unittest
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class ChipseqExecutionContractTests(unittest.TestCase):
    def setUp(self):
        self.contract = json.loads((ROOT / "config/chipseq/execution_contract.json").read_text())
        with (ROOT / "config/chipseq/run_sets.tsv").open(newline="") as handle:
            self.run_sets = list(csv.DictReader(handle, delimiter="\t"))

    def test_chipseq_pin_is_the_released_v120_commit(self):
        self.assertEqual(self.contract["helixforge"], {
            "release": "v1.2.0", "commit": "bd50fe094bc352f698feabe0ae0b8f0e777f5e97",
        })

    def test_peak_paths_complete_without_differential_binding(self):
        composition = self.contract["composition"]
        self.assertEqual(composition["peak_run_modes"], ["consensus", "idr"])
        self.assertIs(composition["postprocess_peak_runs"], True)
        self.assertIs(composition["differential_binding"], False)
        self.assertEqual(composition["peak_run_outputs"], ["annotation", "tracks", "report", "terminal_manifest"])
        for row in self.run_sets:
            if row["run_set"] in {"narrow", "broad"}:
                with self.subTest(study=row["study_accession"], run_set=row["run_set"]):
                    self.assertIn(row["terminal_mode"], composition["peak_run_modes"])

    def test_alignment_only_studies_remain_restricted(self):
        composition = self.contract["composition"]
        self.assertEqual(composition["alignment_only_mode"], "post_alignment")
        self.assertIs(composition["alignment_only_postprocess"], False)
        aligned = {r["study_accession"]: r["terminal_mode"] for r in self.run_sets if r["run_set"] == "alignment"}
        self.assertEqual(aligned, {"PRJNA219440": "post_alignment", "PRJNA384372": "post_alignment", "PRJNA816041": "post_alignment"})

    def test_idr_requires_two_distinct_biological_replicates_and_controls(self):
        for row in self.run_sets:
            if row["terminal_mode"] != "idr":
                continue
            with self.subTest(study=row["study_accession"]):
                self.assertEqual(row["run_set"], "narrow")
                self.assertEqual(row["idr_threshold"], "0.05")
                self.assertEqual(row["idr_rank"], "signal_value")
                path = ROOT / "metadata/chipseq" / row["study_accession"] / row["metadata_file"]
                with path.open(newline="") as handle:
                    records = list(csv.DictReader(handle, delimiter="\t"))
                controls = {r["run_accession"] for r in records if r["is_control"] == "true"}
                groups = defaultdict(list)
                for record in records:
                    if record["is_control"] == "false":
                        self.assertIn(record["control_id"], controls)
                        groups[(record["dataset"], record["condition"], record["target"])].append(record)
                self.assertTrue(groups)
                for records in groups.values():
                    self.assertEqual(len(records), 2)
                    self.assertEqual(len({r["biological_replicate"] for r in records}), 2)

    def test_broad_runs_keep_replicate_support_not_idr(self):
        broad = [r for r in self.run_sets if r["run_set"] == "broad"]
        self.assertEqual(len(broad), 5)
        for row in broad:
            with self.subTest(study=row["study_accession"]):
                self.assertEqual(row["terminal_mode"], "consensus")
                self.assertEqual(row["consensus_policy"], "replicate_support:2")
                self.assertEqual(row["idr_threshold"], "")

    def test_scientific_parameters_are_preserved(self):
        self.assertEqual(self.contract["bam"], {
            "min_mapq": 30, "include_flags": 0, "exclude_flags": 2308,
            "duplicate_mode": "none", "blacklist": None,
        })
        self.assertEqual(self.contract["peaks"], {
            "caller": "macs3", "q_value": 0.01, "duplicate_policy": "all", "mixed_peak_types_split": True,
        })
        self.assertEqual(self.contract["reference"]["genome_size_bp"], 391422299)


if __name__ == "__main__":
    unittest.main()
