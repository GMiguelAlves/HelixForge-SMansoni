from __future__ import annotations

import csv
import hashlib
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HELIXFORGE_V1_COMMIT = "14dc5a75d6c63f20d10c135f0c75138ea76dcc12"
REFERENCE_SHA256 = {
    "genome": "d93a8ff7541d6108b0db088b43e9dc275de45b1d263cd42253877944cceaa1a9",
    "transcriptome": "ef6f9807ba3060d901f3bc89eca6eca2bd0598d162981c7fd3417bc26701dbb4",
    "annotation_gff3": "746fef90d8eb3525d9f2ac0c44f7fd4a8193ac109afe1398deb987a06e6685dd",
    "annotation_gtf": "8b30f5b141ed38b9a99222a18e006880eda38758cd4709a0464141c2a59ccd55",
}
STUDIES = {
    "PRJNA602528": (10, 10),
    "PRJNA597909": (20, 20),
    "PRJEB14695": (138, 23),
    "PRJEB32839": (150, 75),
}


class ProjectContracts(unittest.TestCase):
    def test_registered_study_counts(self) -> None:
        for study, (runs_expected, samples_expected) in STUDIES.items():
            with (ROOT / "metadata" / study / "runs.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                runs = list(csv.DictReader(handle, delimiter="\t"))
            with (ROOT / "metadata" / study / "samples.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                samples = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(runs_expected, len(runs), study)
            self.assertEqual(samples_expected, len(samples), study)

    def test_reference_is_frozen_to_sm_v10_wbps19(self) -> None:
        reference = json.loads(
            (ROOT / "config/reference/reference_manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual("Schistosoma mansoni", reference["organism"]["scientific_name"])
        self.assertEqual("SM_V10", reference["assembly"])
        self.assertEqual("WBPS19", reference["wormbase_parasite_release"])
        self.assertEqual("VERIFIED", reference["status"])
        self.assertEqual("VERIFIED", reference["checksums"]["status"])
        self.assertEqual(
            REFERENCE_SHA256,
            {item["role"]: item["sha256"] for item in reference["artifacts"]},
        )
        self.assertTrue(
            all(item["checksum_status"] == "VERIFIED" for item in reference["artifacts"])
        )

    def test_first_project_is_calibration_not_benchmark(self) -> None:
        state = json.loads(
            (ROOT / "provenance/PRJNA602528/execution_state.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual("FIRST_REAL_SCHISTOSOMA_RNASEQ_PROJECT", state["role"])
        self.assertEqual("OPERATIONAL_CALIBRATION_RUN", state["execution_class"])
        self.assertFalse(state["benchmark"])
        self.assertEqual("NOT_STARTED", state["status"])

    def test_frozen_source_package_checksums(self) -> None:
        frozen = ROOT / "provenance/source_package/frozen"
        entries = []
        for line in (frozen / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines():
            digest, relative = line.split("  ", 1)
            entries.append((digest, relative))
        self.assertEqual(51, len(entries))
        for expected, relative in entries:
            path = frozen / relative
            self.assertTrue(path.is_file(), relative)
            self.assertEqual(expected, hashlib.sha256(path.read_bytes()).hexdigest(), relative)

    def test_no_sequencing_scale_files_are_versioned(self) -> None:
        forbidden = {".bam", ".bai", ".cram", ".crai", ".bigwig", ".bw", ".sif"}
        fastq_suffixes = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
        offenders = []
        for path in ROOT.rglob("*"):
            if not path.is_file():
                continue
            lower = path.name.lower()
            if path.suffix.lower() in forbidden or lower.endswith(fastq_suffixes):
                offenders.append(str(path.relative_to(ROOT)))
        self.assertEqual([], offenders)

    def test_all_execution_states_are_not_started(self) -> None:
        for study in STUDIES:
            state = json.loads(
                (ROOT / "provenance" / study / "execution_state.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual("NOT_STARTED", state["status"], study)
            self.assertEqual("v1.0.0", state["helixforge_release"], study)
            self.assertEqual(HELIXFORGE_V1_COMMIT, state["helixforge_commit"], study)

    def test_helixforge_release_pin_is_exact(self) -> None:
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        server_template = (ROOT / "config/server.env.template").read_text(encoding="utf-8")
        self.assertIn(HELIXFORGE_V1_COMMIT, readme)
        self.assertIn(HELIXFORGE_V1_COMMIT, server_template)

    def test_reference_validation_record_is_consistent(self) -> None:
        validation = json.loads(
            (ROOT / "provenance/reference/reference_validation.json").read_text(
                encoding="utf-8"
            )
        )
        index = json.loads(
            (ROOT / "provenance/reference/salmon_index_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual("PASS", validation["status"])
        self.assertEqual("PASS", validation["contig_compatibility"]["status"])
        self.assertEqual(10960, validation["tx2gene"]["transcriptome_mapped"])
        self.assertEqual(0, validation["tx2gene"]["transcriptome_unmapped"])
        self.assertEqual("VERIFIED_REUSED", index["status"])
        self.assertEqual("1.10.3", index["salmon_version"])
        self.assertEqual(31, index["kmer_size"])

    def test_server_preparation_record_is_sanitized(self) -> None:
        state_path = ROOT / "provenance/server_preparation_state.json"
        state = json.loads(state_path.read_text(encoding="utf-8"))
        self.assertEqual("IMPORT_ONLY", state["prjna602528"]["mode"])
        self.assertEqual("BLOCKED", state["prjna602528"]["differential_expression"])
        self.assertFalse(state["prjna602528"]["fastq_downloaded"])
        self.assertFalse(state["prjna602528"]["scientific_workflow_executed"])

        tracked_text = "\n".join(
            path.read_text(encoding="utf-8", errors="ignore")
            for path in ROOT.rglob("*")
            if path.is_file()
            and ".git" not in path.parts
            and "__pycache__" not in path.parts
            and path.suffix != ".pyc"
        )
        forbidden = (
            "/home/" + "ra236875",
            "/scratch/" + "Schisto-epigenetics",
            "C:\\Users\\" + "dilci",
            "BEGIN OPENSSH " + "PRIVATE KEY",
        )
        for token in forbidden:
            self.assertNotIn(token, tracked_text)

    def test_all_json_documents_parse(self) -> None:
        for path in ROOT.rglob("*.json"):
            json.loads(path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
