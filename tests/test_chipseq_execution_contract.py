from __future__ import annotations

import csv
import hashlib
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

    def test_trimming_is_explicit_and_separate_from_rnaseq(self):
        path = ROOT / "config/chipseq" / self.contract["trimming_policy"]
        policy = json.loads(path.read_text())
        self.assertEqual(policy["status"], "REGISTERED_BEFORE_PRODUCTION")
        self.assertIs(policy["enabled"], True)
        self.assertEqual(policy["tool_version"], "0.6.10")
        self.assertEqual(policy["cutadapt_version"], "5.1")
        self.assertEqual(policy["quality_phred"], 20)
        self.assertEqual(policy["minimum_length_bp"], 25)
        self.assertEqual(policy["adapter_mode"], "auto_detect_per_library")
        self.assertIsNone(policy["adapter1"])
        self.assertIsNone(policy["adapter2"])
        self.assertIs(policy["orphan_rescue"], False)
        self.assertIn("not_exposed", policy["nextseq_poly_g_trimming"])

    def test_bowtie2_reuse_needs_its_own_exact_provenance(self):
        reference = json.loads((ROOT / "config/chipseq/reference_contract.json").read_text())
        self.assertEqual(reference["bowtie2_index_policy"], "require_explicit_manifest_validated_prebuilt_index_no_automatic_rebuild")
        index = reference["bowtie2_index"]
        self.assertEqual(index["aligner"], "bowtie2")
        self.assertEqual(index["bowtie2_version"], "2.5.5")
        self.assertEqual(index["manifest_schema_version"], "1.2")
        self.assertEqual(index["basename"], "genome")
        self.assertEqual(index["extra_args"], [])
        self.assertEqual(index["location_policy"], "private_site_configuration_only")
        self.assertEqual(index["validation_state"], "PASS_SLURM_BUILD_AND_FRESH_REUSE")
        self.assertEqual(reference["fasta_sha256"], "64061ee54997f652dccf1a05718ae941843b762aed4bf96b0c0f699ecbf01fb3")
        self.assertEqual(reference["gtf_sha256"], "0ba2e3ff9a6328080a20f7b91bc7aa66c558eda5fdb6c382349d5e1f1a5f2d0e")

    def test_index_validation_record_matches_the_frozen_bundle(self):
        reference = json.loads((ROOT / "config/chipseq/reference_contract.json").read_text())
        index = reference["bowtie2_index"]
        record = json.loads((ROOT / "config/chipseq" / index["validation_record"]).read_text())
        self.assertEqual(record["status"], "PASS")
        self.assertEqual(record["helixforge_commit"], self.contract["helixforge"]["commit"])
        self.assertEqual(record["reference_sha256"], reference["fasta_sha256"])
        self.assertEqual(record["annotation_sha256"], reference["gtf_sha256"])
        self.assertEqual(record["bowtie2_version"], index["bowtie2_version"])
        self.assertEqual(record["fresh_reuse"], {"status": "PASS", "resume_used": False, "build_tasks": 0, "validation_tasks": 1})
        self.assertEqual(record["construction"], {"status": "PASS", "build_tasks": 1, "validation_tasks": 0})
        self.assertEqual(len(record["files"]), 6)
        digest = hashlib.sha256(json.dumps(record["files"], sort_keys=True).encode()).hexdigest()
        self.assertEqual(digest, index["bundle_sha256"])
        self.assertEqual(digest, record["bundle_sha256"])


if __name__ == "__main__":
    unittest.main()
