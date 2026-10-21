"""Native-source conservation and fail-closed admission checks for atlas ChIP."""
import copy
import csv
import json
import gzip
import io
import hashlib
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts/atlas"))
from hf_chipseq_evidence import load_hf_chipseq
from hf_chipseq_regional import load_regional_study
from build_rnaseq_atlas import version_report_assets


class NativeChipTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        (ROOT / ".codex_tmp").mkdir(exist_ok=True)
        atlas = json.loads((ROOT / "analyses/atlas/atlas_config.json").read_text(encoding="utf-8"))
        cls.config = {**atlas["chipseq_integration"], "reference_id": atlas["reference_id"], "studies": [atlas["chipseq_integration"]["studies"][0]]}
        with (ROOT / atlas["gene_annotations"]).open(encoding="utf-8") as handle:
            cls.annotations = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

    def test_all_native_segments_and_gene_assignments_are_preserved(self):
        result = load_hf_chipseq(ROOT, self.config, self.annotations)
        self.assertEqual(len(result["segments"]), 2683)
        self.assertEqual(sum(r["category"] == "intergenic" for r in result["segments"]), 1458)
        self.assertEqual(len(result["loci"]), 2683)
        self.assertEqual(len({r["gene_id"] for r in result["loci"]}), 2014)
        self.assertEqual(sum(r["n_peaks"] for r in result["evidence"]), 2683)
        self.assertTrue(all(r["support"] >= 1 for r in result["segments"]))
        self.assertEqual(sum(r["n_samples"] for r in result["metadata"]), 6)
        self.assertEqual(sum(r["n_controls"] for r in result["metadata"]), 2)
        native = {r["peak_id"]: r for r in result["segments"]}
        for row in result["loci"]:
            segment = native[row["peak_id"]]
            self.assertEqual((row["peak_start"], row["peak_end"]), (int(segment["start"]), int(segment["end"])))
            self.assertEqual(row["peak_location"], segment["category"])
        profile = result["profiles"][0]
        self.assertEqual(profile["qc"]["differential_binding"], "NOT_ELIGIBLE_IN_CURRENT_ANALYSIS")
        self.assertTrue(all(p["provider"] == "homer_4_11" for p in profile["annotation_policies"]))
        self.assertEqual(sum(r["support"] == 1 for r in result["segments"]),2382)
        self.assertEqual(len(result["core_segments"]),465)
        self.assertEqual(len(result["regional_support"]),8049)
        self.assertTrue(all(r["gene_ids"] for r in result["segments"] if r["category"] == "intergenic"))
        self.assertTrue(all(r["gene_assignment_ambiguous"] == "not_assessed" for r in result["loci"]))

    def test_incompatible_reference_or_missing_genes_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Incompatible"):
            load_hf_chipseq(ROOT, {**self.config, "reference_id": "another_reference"}, self.annotations)
        missing = dict(self.annotations)
        missing.pop("Smp_244880")
        with self.assertRaisesRegex(ValueError, "outside the (atlas|reference)"):
            load_hf_chipseq(ROOT, self.config, missing)

    def test_changed_source_cannot_silently_enter_the_atlas(self):
        entry = self.config["studies"][0]
        source = ROOT / entry["source_manifest"]
        with tempfile.TemporaryDirectory(dir=ROOT / ".codex_tmp") as directory:
            root = Path(directory)
            self.assertTrue(root.resolve().is_relative_to((ROOT / ".codex_tmp").resolve()))
            copied = root / entry["source_manifest"]
            shutil.copytree(source.parent, copied.parent)
            native = copied.parent / "results/chipseq/PRJNA1090249/broad/chipseq/peak_annotation/cercaria/peak_gene_associations.tsv.gz"
            native.write_bytes(native.read_bytes()+b"\n")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                load_hf_chipseq(root, self.config, self.annotations)
        changed = copy.deepcopy(self.config)
        changed["studies"][0]["commit"] = "0" * 40
        with self.assertRaisesRegex(ValueError, "frozen study/commit"):
            load_hf_chipseq(ROOT, changed, self.annotations)

    def test_legacy_atomic_snapshot_remains_distinct(self):
        config = copy.deepcopy(self.config)
        config["studies"][0].update(commit="a46707208b6ac326f6836a77c353f86319e5f550",source_manifest="analyses/atlas/local_inputs/hf_chipseq/PRJNA1090249/source_manifest.json")
        result=load_hf_chipseq(ROOT,config,self.annotations)
        self.assertEqual(len(result["segments"]),465)
        self.assertEqual(len(result["core_segments"]),0)
        self.assertEqual(len(result["library_peaks"]),3080)

    def mutate_regional_table(self,role,change):
        entry=self.config["studies"][0]
        source_path=ROOT/entry["source_manifest"]
        with tempfile.TemporaryDirectory(dir=ROOT/".codex_tmp") as folder:
            base=Path(folder)/"snapshot"
            shutil.copytree(source_path.parent,base)
            result=base/"results/chipseq/PRJNA1090249"
            manifest_path=result/"broad/chipseq/consensus/cercaria/manifest.json"
            manifest=json.loads(manifest_path.read_text())
            artifact=manifest["artifacts"][role]
            path=result/artifact["path"]
            with gzip.open(path,"rt",encoding="utf-8") as handle:
                reader=csv.DictReader(handle,delimiter="\t");fields=reader.fieldnames;rows=list(reader)
            rows=change(rows)
            buffer=io.StringIO();writer=csv.DictWriter(buffer,fieldnames=fields,delimiter="\t",lineterminator="\n");writer.writeheader();writer.writerows(rows)
            path.write_bytes(gzip.compress(buffer.getvalue().encode()))
            artifact["sha256"]=hashlib.sha256(path.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest))
            source=json.loads(source_path.read_text())
            return load_regional_study(ROOT,result,entry,source,self.config,self.annotations)

    def test_zero_support_library_rows_cannot_be_dropped(self):
        def drop(rows):
            index=next(i for i,r in enumerate(rows) if r["called_peak_count"]=="0")
            return rows[:index]+rows[index+1:]
        with self.assertRaisesRegex(ValueError,"omits libraries"):
            self.mutate_regional_table("support_by_library.tsv",drop)

    def test_unassessed_ambiguity_cannot_become_false(self):
        def change(rows):
            rows[0]["gene_assignment_ambiguous"]="false";return rows
        with self.assertRaisesRegex(ValueError,"ambiguity state"):
            self.mutate_regional_table("annotated_peaks.tsv",change)

    def test_changed_scientific_figure_gets_a_new_browser_url(self):
        with tempfile.TemporaryDirectory(dir=ROOT / ".codex_tmp") as directory:
            output = Path(directory).resolve()
            self.assertTrue(output.is_relative_to((ROOT / ".codex_tmp").resolve()))
            (output / "figures").mkdir()
            figure = output / "figures/study.svg"
            figure.write_text('<svg>historical</svg>', encoding="utf-8")
            before = version_report_assets('<img src="figures/study.svg">', output)
            figure.write_text('<svg>native HP1 results</svg>', encoding="utf-8")
            after = version_report_assets('<img src="figures/study.svg">', output)
            self.assertNotEqual(before, after)
            self.assertIn('.svg?v=', after)
            self.assertEqual(version_report_assets('<img src="figures/${study}.svg">', output), '<img src="figures/${study}.svg">')


if __name__ == "__main__":
    unittest.main()
