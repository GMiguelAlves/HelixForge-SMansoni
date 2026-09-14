from __future__ import annotations

import csv
import hashlib
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HELIXFORGE_V1_COMMIT = "14dc5a75d6c63f20d10c135f0c75138ea76dcc12"
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
        self.assertEqual("PENDING_VERIFIED_DOWNLOAD", reference["checksums"]["status"])

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


if __name__ == "__main__":
    unittest.main()
