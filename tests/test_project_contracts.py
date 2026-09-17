from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HELIXFORGE_V1_COMMIT = "e41d221657b8e0bf2700bccd547e15032ccac36f"
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


def load_prjeb14695_preflight():
    path = ROOT / "scripts/validate/preflight_prjeb14695.py"
    spec = importlib.util.spec_from_file_location("preflight_prjeb14695", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def write_ena_fixture(root: Path, destination: Path) -> None:
    runs = read_tsv(root / "metadata/PRJEB14695/runs.tsv")
    samples = {row["sample_id"]: row for row in read_tsv(root / "metadata/PRJEB14695/samples.tsv")}
    fields = [
        "run_accession", "study_accession", "secondary_study_accession",
        "experiment_accession", "sample_accession", "secondary_sample_accession",
        "sample_alias", "library_name", "library_layout", "library_strategy",
        "library_source", "fastq_ftp", "fastq_md5", "fastq_bytes",
        "read_count", "base_count", "first_public",
    ]
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for run in runs:
            sample = samples[run["sample_id"]]
            writer.writerow({
                "run_accession": run["run_accession"],
                "study_accession": "PRJEB14695",
                "secondary_study_accession": "ERP016356",
                "experiment_accession": run["experiment_accession"],
                "sample_accession": run["biosample"],
                "secondary_sample_accession": run["secondary_biosample"],
                "sample_alias": sample["source_sample_name"],
                "library_name": run["library_name"],
                "library_layout": "PAIRED",
                "library_strategy": "RNA-Seq",
                "library_source": "TRANSCRIPTOMIC",
                "fastq_ftp": ";".join((run["fastq_1_url"], run["fastq_2_url"])),
                "fastq_md5": ";".join((run["fastq_1_md5"], run["fastq_2_md5"])),
                "fastq_bytes": ";".join((run["fastq_1_bytes"], run["fastq_2_bytes"])),
                "read_count": run["read_count"],
                "base_count": run["base_count"],
                "first_public": run["ena_first_public"],
            })


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ProjectContracts(unittest.TestCase):
    def test_prjeb14695_scientific_preflight_materializes_23_observations(self) -> None:
        preflight = load_prjeb14695_preflight()
        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            ena = temporary_root / "ena.tsv"
            output = temporary_root / "output"
            write_ena_fixture(ROOT, ena)
            report = preflight.validate(ROOT, ena, output)
            self.assertEqual("PASS", report["status"])
            self.assertEqual("PASS", report["gates"]["run_to_sample_mapping"])
            self.assertEqual("PASS", report["gates"]["technical_run_aggregation"])
            self.assertEqual(138, report["run_inventory"]["mapped"])
            self.assertEqual(23, report["sample_accounting"]["deseq2_observations"])
            self.assertEqual(6, report["sample_accounting"]["runs_per_sample_min"])
            self.assertEqual(6, report["sample_accounting"]["runs_per_sample_max"])
            self.assertTrue(report["design"]["full_rank"])
            self.assertTrue(report["design"]["all_contrasts_estimable"])
            self.assertEqual(23, len(read_tsv(output / "technical_run_aggregation.tsv")))
            self.assertEqual(138, len(read_tsv(output / "run_to_sample_mapping.tsv")))

    def test_prjeb14695_scientific_preflight_rejects_conflicting_assignment(self) -> None:
        preflight = load_prjeb14695_preflight()
        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            project_root = temporary_root / "project"
            shutil.copytree(ROOT / "metadata/PRJEB14695", project_root / "metadata/PRJEB14695")
            shutil.copytree(ROOT / "config/PRJEB14695", project_root / "config/PRJEB14695")
            ena = temporary_root / "ena.tsv"
            output = temporary_root / "output"
            write_ena_fixture(ROOT, ena)
            metadata_path = project_root / "metadata/PRJEB14695/metadata.csv"
            with metadata_path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
                fields = list(rows[0])
            conflicting = next(row for row in rows if row["sample_id"] != rows[0]["sample_id"])
            rows[0]["sample_id"] = conflicting["sample_id"]
            with metadata_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            report = preflight.validate(project_root, ena, output)
            self.assertEqual("FAIL", report["status"])
            self.assertEqual("FAIL", report["gates"]["run_to_sample_mapping"])
            self.assertTrue(any("run-to-sample conflicts" in error for error in report["errors"]))

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
        self.assertIn(
            state["status"],
            {
                "NOT_STARTED",
                "ACQUISITION_IN_PROGRESS",
                "BLOCKED_PRE_WORKFLOW_INDEX_REUSE_CONFLICT",
                "READY_FOR_WORKFLOW",
                "WORKFLOW_RUNNING",
                "WORKFLOW_COMPLETE",
                "READY_FOR_REVIEW",
                "ACCEPTED",
            },
        )
        if state["status"] != "NOT_STARTED":
            self.assertEqual("PASS", state["acquisition"]["integrity"])

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

    def test_unstarted_projects_remain_frozen(self) -> None:
        for study in set(STUDIES) - {"PRJNA602528", "PRJNA597909"}:
            state = json.loads(
                (ROOT / "provenance" / study / "execution_state.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual("NOT_STARTED", state["status"], study)
            self.assertEqual("v1.0.1", state["helixforge_release"], study)
            self.assertEqual(HELIXFORGE_V1_COMMIT, state["helixforge_commit"], study)

    def test_prjna597909_execution_is_complete(self) -> None:
        state = json.loads(
            (ROOT / "provenance/PRJNA597909/execution_state.json").read_text(
                encoding="utf-8"
            )
        )
        preflight = json.loads(
            (ROOT / "provenance/PRJNA597909/preflight.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual("READY_FOR_REVIEW", state["status"])
        self.assertEqual("COMPLETE", state["phase"])
        self.assertTrue(state["design"]["full_rank"])
        self.assertTrue(state["design"]["contrasts_estimable"])
        self.assertEqual("PASS", state["acquisition"]["integrity"])
        self.assertEqual("PASS_WITH_LIMITATIONS", state["execution"]["process_status"])
        self.assertEqual("TERMINAL_MANIFEST_ONLY", state["execution"]["reentry_scope"])
        self.assertFalse(
            state["execution"]["scientific_stages_recomputed_during_recovery"]
        )
        self.assertEqual("PASS", state["storage"]["cleanup"])
        self.assertEqual(0, state["storage"]["scratch_bytes_after_cleanup"])
        self.assertEqual("PASS_WITH_LIMITATIONS", state["validation"]["analysis"])
        self.assertFalse(state["safety"]["next_project_authorized"])
        self.assertEqual("PASS", preflight["status"])
        self.assertEqual("PASS", preflight["gates"]["salmon_index_reuse"])
        self.assertEqual("PASS", preflight["gates"]["report_contract_preflight"])

    def test_launcher_requires_validated_prebuilt_index(self) -> None:
        launcher = (ROOT / "scripts/run_study.sh").read_text(encoding="utf-8")
        self.assertIn("--salmon_prebuilt_index", launcher)
        self.assertIn("--salmon_prebuilt_index_manifest", launcher)
        self.assertNotIn("salmon index ", launcher)

    def test_launcher_uses_resume_only_when_explicitly_requested(self) -> None:
        launcher = (ROOT / "scripts/run_study.sh").read_text(encoding="utf-8")
        self.assertIn('"${HF_RESUME:-0}" == 1', launcher)
        self.assertIn("args+=(-resume)", launcher)
        self.assertNotIn("-w \"$HF_WORK_ROOT\" -resume", launcher)

    def test_launcher_requires_certified_python_runtime(self) -> None:
        launcher = (ROOT / "scripts/run_study.sh").read_text(encoding="utf-8")
        preflight = (ROOT / "scripts/validate/preflight_server.sh").read_text(
            encoding="utf-8"
        )
        slurm = (ROOT / "config/slurm.config").read_text(encoding="utf-8")
        template = (ROOT / "config/server.env.template").read_text(encoding="utf-8")
        for source in (launcher, preflight, template):
            self.assertIn("HF_CERTIFIED_RUNTIME_PATH", source)
        self.assertIn("import jsonschema, sys", launcher)
        self.assertIn("import jsonschema, sys", preflight)
        self.assertIn('env PATH="$HF_CERTIFIED_RUNTIME_PATH"', launcher)
        self.assertIn("env.PATH = hfEnv.HF_CERTIFIED_RUNTIME_PATH", slurm)
        self.assertNotIn("--export=ALL", slurm)

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

        tracked = subprocess.run(
            ["git", "-c", f"safe.directory={ROOT.as_posix()}", "-C", str(ROOT), "ls-files", "-z"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.split("\0")
        tracked_text = "\n".join(
            (ROOT / relative).read_text(encoding="utf-8", errors="ignore")
            for relative in tracked
            if relative and (ROOT / relative).is_file()
        )
        forbidden = (
            "/home/" + "ra" + "236875",
            "/scratch/" + "Schisto-epigenetics",
            "C:\\Users\\" + "dilci",
            "BEGIN OPENSSH " + "PRIVATE KEY",
        )
        for token in forbidden:
            self.assertNotIn(token, tracked_text)

    def test_prjna602528_import_only_results_are_complete(self) -> None:
        result_root = ROOT / "results/PRJNA602528"
        manifest_path = result_root / "manifests/rnaseq_run_manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual("rnaseq_run_manifest", manifest["type"])
        self.assertEqual("complete", manifest["status"])
        self.assertEqual("salmon", manifest["quantification_method"])
        self.assertEqual([], manifest["contrasts"])
        self.assertEqual(10, len(manifest["samples"]))
        self.assertEqual(12, len(manifest["artifacts"]))

        for artifact in manifest["artifacts"]:
            self.assertEqual("manifest_relative", artifact["location"]["kind"])
            path = (manifest_path.parent / artifact["location"]["path"]).resolve()
            self.assertTrue(path.is_file(), artifact["artifact_id"])
            observed = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(artifact["checksum"]["value"], observed)

        with (result_root / "expression/counts_matrix.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        self.assertEqual(10, len(rows[0]) - 1)
        self.assertEqual(9914, len(rows) - 1)

        qc_rows = read_tsv(result_root / "qc/PRJNA602528_qc_summary.tsv")
        self.assertEqual(10, len(qc_rows))
        self.assertEqual(1, sum(row["qc_flag"] == "REVIEW" for row in qc_rows))

        report = (result_root / "reports/PRJNA602528_execution_report.md").read_text(
            encoding="utf-8"
        )
        self.assertIn("DIFFERENTIAL_EXPRESSION = NOT_PERFORMED_BY_DESIGN", report)
        self.assertTrue(report.rstrip().endswith("READY_FOR_PRJNA602528_REVIEW"))

        state = json.loads(
            (ROOT / "provenance/PRJNA602528/execution_state.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual("NOT_APPLICABLE", state["validation"]["differential_expression"])
        self.assertEqual("PASS", state["validation"]["cleanup"])
        self.assertEqual("PASS_WITH_LIMITATIONS", state["validation"]["analysis"])

    def test_prjna597909_full_results_are_complete(self) -> None:
        result_root = ROOT / "results/PRJNA597909"
        manifest_path = result_root / "rnaseq/rnaseq_run_manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual("rnaseq_run_manifest", manifest["type"])
        self.assertEqual("complete", manifest["status"])
        self.assertEqual("salmon", manifest["quantification_method"])
        self.assertEqual(20, len(manifest["samples"]))
        self.assertEqual(4, len(manifest["contrasts"]))

        portable_artifacts = [
            artifact
            for artifact in manifest["artifacts"]
            if artifact["location"]["kind"] == "manifest_relative"
        ]
        self.assertEqual(8, len(portable_artifacts))
        for artifact in portable_artifacts:
            path = (manifest_path.parent / artifact["location"]["path"]).resolve()
            self.assertTrue(path.is_file(), artifact["artifact_id"])
            observed = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(artifact["checksum"]["value"], observed)

        validation = json.loads(
            (result_root / "manifests/terminal_manifest_validation.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual("complete", validation["status"])
        self.assertEqual("valid", validation["schema"])
        self.assertEqual("valid", validation["semantic"])

        qc_rows = read_tsv(result_root / "qc/PRJNA597909_qc_summary.tsv")
        self.assertEqual(20, len(qc_rows))
        self.assertTrue(all(row["classification"] == "PASS" for row in qc_rows))

        for matrix in ("counts_matrix.tsv", "tpm_matrix.tsv", "length_matrix.tsv"):
            with (result_root / "expression" / matrix).open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.reader(handle, delimiter="\t"))
            self.assertEqual(20, len(rows[0]) - 1, matrix)
            self.assertEqual(9914, len(rows) - 1, matrix)

        expected_de = {
            "condition__juvenile_pzq_vs_juvenile_control": (204, 107, 97),
            "condition__adult_pzq_vs_adult_control": (636, 338, 298),
            "condition__adult_control_vs_juvenile_control": (1131, 757, 374),
            "condition__adult_pzq_vs_juvenile_pzq": (1075, 763, 312),
        }
        contrasts = read_tsv(
            result_root / "differential_expression/contrast_summary.tsv"
        )
        self.assertEqual(4, len(contrasts))
        for row in contrasts:
            self.assertEqual(
                expected_de[row["contrast"]],
                (int(row["significant"]), int(row["up"]), int(row["down"])),
            )
            self.assertEqual("PASS", row["status"])

        report = (result_root / "reports/PRJNA597909_execution_report.md").read_text(
            encoding="utf-8"
        )
        self.assertIn("PRJNA597909_RNASEQ_ANALYSIS = PASS_WITH_LIMITATIONS", report)
        self.assertTrue(report.rstrip().endswith("READY_FOR_PRJNA597909_REVIEW"))

        final = json.loads(
            (result_root / "reports/final_validation.json").read_text(encoding="utf-8")
        )
        self.assertEqual("PASS_WITH_LIMITATIONS", final["status"])
        self.assertTrue(final["ready_for_review"])
        self.assertFalse(final["next_project_authorized"])

        expected_html = {
            "PRJNA597909_multiqc.html",
            "nextflow_dag.html",
            "nextflow_execution_report.html",
            "nextflow_timeline.html",
        }
        observed_html = {
            path.name for path in (result_root / "reports").glob("*.html")
        }
        self.assertEqual(expected_html, observed_html)
        self.assertEqual([], list(result_root.rglob("*_fastqc.html")))

    def test_all_json_documents_parse(self) -> None:
        for path in ROOT.rglob("*.json"):
            json.loads(path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
