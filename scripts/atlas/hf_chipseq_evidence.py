"""Import reviewed HelixForge ChIP results, preserving study and source classes.

No peak calling, reassignment or differential binding is performed here.
The CLI freezes small published artifacts from a specific Git commit.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from hf_chipseq_regional import load_regional_study, SOURCE_FIELDS
from hf_chipseq_catalog import load_catalog


REPOSITORY = "https://github.com/GMiguelAlves/HelixForge-SMansoni"
METADATA_FIELDS = ["study_id", "context_id", "mark_or_factor", "stage_or_condition", "n_samples", "n_controls", "replicates", "batches", "qc_status"]
EVIDENCE_FIELDS = ["study_id", "context_id", "gene_id", "functional_name", "mark_or_factor", "stage_or_condition", "n_peaks", "n_promoter_peaks", "peak_locations", "regulatory_class", "expected_effect"]
LOCUS_FIELDS = ["study_id", "context_id", "gene_id", "stage_or_condition", "mark_or_factor", "peak_id", "peak_chrom", "peak_start", "peak_end", "peak_midpoint", "peak_location", "promoter_flag", "position_class", "support"]
POSITION_FIELDS = ["study_id", "context_id", "gene_id", "stage_or_condition", "mark_or_factor", "position_class", "n_unique_peaks"]
METADATA_FIELDS += ["consensus_method", "annotation_provider", "min_libraries"]
EVIDENCE_FIELDS += ["n_single_library", "n_multiple_libraries", "n_ambiguous", "n_ambiguity_not_assessed"]
LOCUS_FIELDS += SOURCE_FIELDS

LIMITATIONS = (
    "Descriptive HP1 occupancy only. Stage and strain are confounded "
    "(cercaria SmVEN; sporocyst SmDFO), and independence of the three deposited "
    "IP libraries per group is not established. FRiP is 0.201–3.161%. The 200 bp "
    "model-free extension is an operational assumption, not a measured fragment "
    "length. Consensus intervals are support segments (at least two deposited "
    "IP libraries), not verified broad domains. Differential binding was not run."
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def freeze_source(repository: Path, revision: str, study: str, output: Path) -> None:
    commit = subprocess.check_output(
        ["git", "rev-parse", f"{revision}^{{commit}}"], cwd=repository, text=True
    ).strip()
    prefix = f"results/chipseq/{study}/"
    files = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", commit], cwd=repository, text=True
    ).splitlines()
    exact = {
        f"{prefix}README.md", f"{prefix}qc/qc_summary.json",
        f"{prefix}qc/qc_sample_summary.tsv", f"{prefix}manifests/publication_validation.json",
        f"{prefix}manifests/reference_bundle.manifest.json",
        f"{prefix}reports/{study}_execution_report.md",
        f"{prefix}broad/chipseq/chipseq_run_manifest.json",
        f"{prefix}broad/chipseq/consensus/consolidation_summary.tsv",
        f"{prefix}broad/chipseq/consensus/consolidation_summary.json",
        f"{prefix}broad/chipseq/peak_qc/peak_qc_summary.tsv",
        f"provenance/chipseq/{study}/execution_provenance.json",
        "docs/chipseq/prjna1090249_peak_policy.md",
        "docs/chipseq/differential_binding_review.md",
    }
    selected = [name for name in files if name in exact or (
        name.startswith(prefix + "broad/chipseq/peak_annotation/")
        and (name.endswith(".tsv") or name.endswith("/manifest.json"))
    ) or (
        name.startswith(prefix + "broad/chipseq/consensus/")
        and Path(name).name in {"consolidated_peaks.tsv", "consolidated_peaks.bed", "manifest.json", "statistics.json", "replicate_evidence.tsv"}
    )]
    regional = f"{prefix}manifests/regional_publication.json" in files
    if regional:
        exact -= {f"{prefix}broad/chipseq/chipseq_run_manifest.json", f"{prefix}broad/chipseq/consensus/consolidation_summary.tsv",
                  f"{prefix}broad/chipseq/consensus/consolidation_summary.json", f"{prefix}broad/chipseq/peak_qc/peak_qc_summary.tsv"}
        exact |= {f"{prefix}manifests/regional_publication.json", f"{prefix}broad/regional_summary.tsv"}
        selected = [name for name in files if name in exact or
                    name.startswith(prefix + "broad/chipseq/consensus/") and name.endswith((".gz", "/manifest.json", "/statistics.json")) or
                    name.startswith(prefix + "broad/chipseq/peak_annotation/") and name.endswith(".gz") or
                    name.startswith(prefix + "broad/080-peak-calling/") and (name.endswith("/peaks.broadPeak") or name.endswith("/manifest.json"))]
    if not exact.issubset(files):
        raise ValueError(f"Missing published artifacts: {sorted(exact - set(files))}")
    records = []
    for name in selected:
        content = subprocess.check_output(["git", "show", f"{commit}:{name}"], cwd=repository)
        target = output / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
        records.append({"path": name, "sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)})
    manifest = {"schema_version": "1.0", "study_id": study, "repository": REPOSITORY,
                "commit": commit, "scope": "small_published_HF_artifacts_no_reads_or_tracks", "representation": "nfcore_2_1_min1_homer_4_11" if regional else "atomic_support", "files": records}
    (output / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"commit": commit, "files": len(records), "bytes": sum(r["bytes"] for r in records)}, indent=2))


def load_hf_chipseq(root: Path, config: dict, annotations: dict) -> dict:
    metadata, evidence, positions, loci, segments, profiles, inputs, library_peaks = [], [], [], [], [], [], [], []
    core_segments, regional_support, idr_segments = [], [], []
    ids = [entry["id"] for entry in config["studies"]]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate ChIP study in the atlas configuration")
    for entry in config["studies"]:
        source_manifest = root / entry["source_manifest"]
        source = json.loads(source_manifest.read_text(encoding="utf-8"))
        if source["commit"] != entry["commit"] or source["study_id"] != entry["id"]:
            raise ValueError("ChIP source does not match the configured frozen study/commit")
        inputs.append(source_manifest)
        base = source_manifest.parent
        for record in source["files"]:
            path = base / record["path"]
            if not path.is_file() or digest(path) != record["sha256"]:
                raise ValueError(f"ChIP source checksum mismatch: {path}")
            inputs.append(path)
        study = entry["id"]
        result = base / f"results/chipseq/{study}"
        if entry.get("admission_scope") == "published_descriptive_catalog_only":
            catalog = load_catalog(root, result, entry, source, config, annotations)
            for target, key in [(metadata,"metadata"),(evidence,"evidence"),(positions,"positions"),(loci,"loci"),(segments,"segments"),(profiles,"profiles"),(library_peaks,"library_peaks"),(core_segments,"core_segments"),(regional_support,"regional_support"),(idr_segments,"idr_segments")]:
                target.extend(catalog[key])
            continue
        if (result / "manifests/regional_publication.json").is_file():
            regional = load_regional_study(root, result, entry, source, config, annotations)
            for target, key in [(metadata,"metadata"),(evidence,"evidence"),(positions,"positions"),(loci,"loci"),(segments,"segments"),(profiles,"profiles"),(library_peaks,"library_peaks"),(core_segments,"core_segments"),(regional_support,"regional_support")]:
                target.extend(regional[key])
            continue
        chip = result / "broad/chipseq"
        qc = json.loads((result / "qc/qc_summary.json").read_text(encoding="utf-8"))
        if (qc["atlas_inclusion"], qc["technical_execution"], qc["processing_status"]) != ("ACCEPTED", "PASS", "COMPLETE"):
            raise ValueError(f"ChIP study has not been accepted: {study}")
        samples = read_tsv(result / "qc/qc_sample_summary.tsv")
        groups = read_tsv(chip / "consensus/consolidation_summary.tsv")
        annotated = read_tsv(chip / "peak_annotation/peak_annotation_aggregate/annotated_peaks.tsv")
        associations = read_tsv(chip / "peak_annotation/peak_annotation_aggregate/peak_gene_associations.tsv")
        peak_lookup = {row["peak_id"]: row for row in annotated}
        if len(peak_lookup) != len(annotated):
            raise ValueError("Duplicate consensus segment IDs in native annotation")
        native_segments, contexts, policies = {}, {}, []
        for group in groups:
            if group["genome_id"] != config["reference_id"] or group["status"] != "complete":
                raise ValueError("Incompatible/incomplete native ChIP consensus")
            condition, mark = group["condition"], group["target"]
            context = f"{study}|{condition}|{mark}"
            contexts[group["group_id"]] = (condition, mark, context)
            matching = [r for r in samples if r["condition"] == condition]
            ips = [r for r in matching if r["is_control"].lower() == "false"]
            controls = [r for r in matching if r["is_control"].lower() == "true"]
            metadata.append({"study_id": study, "context_id": context, "mark_or_factor": mark,
                             "stage_or_condition": condition, "n_samples": len(ips), "n_controls": len(controls),
                             "replicates": ";".join(r["sample_id"] for r in ips),
                             "batches": "not_established", "qc_status": qc["qc_classification"]})
            consensus = chip / "consensus" / group["group_id"] / (group["group_id"] + ".replicate_support.consensus_result")
            replicate_source = consensus / "replicate_evidence.tsv"
            if not replicate_source.is_file():
                raise ValueError("Frozen ChIP source lacks per-library peaks; refresh the source snapshot")
            replicate_manifest = json.loads((consensus / "manifest.json").read_text(encoding="utf-8"))
            replicates = {r["evidence_replicate_id"]: r for r in replicate_manifest["replicates"]}
            observed_counts = Counter()
            for peak in read_tsv(replicate_source):
                rep = replicates[peak["replicate_id"]]
                if not 0 <= int(peak["start"]) < int(peak["end"]):
                    raise ValueError("Invalid per-library ChIP peak coordinates")
                observed_counts[rep["sample_id"]] += 1
                library_peaks.append({"study_id": study, "stage_or_condition": condition,
                    "mark_or_factor": mark, "sample_id": rep["sample_id"], "run_accession": rep["record_id"],
                    "peak_id": peak["original_peak_name"], "chrom": peak["chrom"],
                    "start": int(peak["start"]), "end": int(peak["end"]),
                    "signal_value": float(peak["signal_value"]), "p_value": float(peak["p_value"]),
                    "q_value": float(peak["q_value"])})
            if any(observed_counts[r["sample_id"]] != int(r["peak_count"]) for r in ips):
                raise ValueError("Per-library peaks differ from source QC counts")
            for row in read_tsv(consensus / "consolidated_peaks.tsv"):
                if int(row["support"]) < 2 or int(row["end"]) <= int(row["start"]) or row["peak_id"] in native_segments:
                    raise ValueError("Invalid consensus segment/support")
                native_segments[row["peak_id"]] = row
            policy_path = chip / "peak_annotation" / (group["group_id"] + ".annotation.peak_annotation") / "manifest.json"
            policy = json.loads(policy_path.read_text(encoding="utf-8"))
            if policy["build"] != config["reference_id"]:
                raise ValueError("ChIP annotation build differs from atlas")
            policies.append({"condition": condition, **policy["parameters"]})
        if native_segments.keys() != peak_lookup.keys():
            raise ValueError("Native annotation does not conserve consensus segment identities")
        for row in annotated:
            native = native_segments[row["peak_id"]]
            if any(str(row[k]) != str(native[k]) for k in ("chrom", "start", "end")):
                raise ValueError("Native consensus/annotation coordinate mismatch")
            stage, mark, context = contexts[row["source_id"]]
            segments.append({"study_id": study, "context_id": context, "stage_or_condition": stage,
                             "mark_or_factor": mark, **{k: row[k] for k in ("peak_id", "chrom", "start", "end", "category", "gene_ids")},
                             "support": int(native["support"]), "support_libraries": native["support_replicates"]})
        summary, position_sets, seen = {}, defaultdict(set), set()
        for row in associations:
            gene, peak = row["gene_id"], peak_lookup[row["peak_id"]]
            if gene not in annotations:
                raise ValueError(f"Native ChIP gene is outside the atlas: {gene}")
            stage, mark, context = contexts[row["source_id"]]
            key = (gene, context, row["peak_id"])
            if key in seen or peak["category"] != row["category"] or annotations[gene]["chromosome"] != peak["chrom"]:
                raise ValueError("Duplicate/inconsistent native gene–segment association")
            seen.add(key)
            category = row["category"]
            position = "promoter_TSS" if category == "promoter" else "gene_body" if category in {"gene", "exon", "intron"} else "other_annotated"
            record = summary.setdefault((gene, context), {"study_id": study, "context_id": context,
                "gene_id": gene, "functional_name": annotations[gene]["functional_name"], "mark_or_factor": mark,
                "stage_or_condition": stage, "n_peaks": 0, "n_promoter_peaks": 0, "peak_locations": set(),
                "regulatory_class": "observed_factor_association", "expected_effect": "not_inferred_from_occupancy"})
            record["n_peaks"] += 1
            record["n_promoter_peaks"] += int(category == "promoter")
            record["peak_locations"].add(category)
            position_sets[(gene, context, stage, mark, position)].add(row["peak_id"])
            loci.append({"study_id": study, "context_id": context, "gene_id": gene,
                         "stage_or_condition": stage, "mark_or_factor": mark, "peak_id": row["peak_id"],
                         "peak_chrom": peak["chrom"], "peak_start": int(peak["start"]), "peak_end": int(peak["end"]),
                         "peak_midpoint": (int(peak["start"]) + int(peak["end"])) // 2,
                         "peak_location": category, "promoter_flag": str(category == "promoter").lower(),
                         "position_class": position, "support": int(native_segments[row["peak_id"]]["support"])})
        for row in summary.values():
            row["peak_locations"] = ";".join(sorted(row["peak_locations"]))
            evidence.append(row)
        for (gene, context, stage, mark, position), ids in position_sets.items():
            positions.append({"study_id": study, "context_id": context, "gene_id": gene,
                              "stage_or_condition": stage, "mark_or_factor": mark,
                              "position_class": position, "n_unique_peaks": len(ids)})
        provenance_path = base / f"provenance/chipseq/{study}/execution_provenance.json"
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        profile = {"study_id": study, "reference_id": config["reference_id"], "target": groups[0]["target"],
                   "repository": source["repository"], "commit": source["commit"], "qc": qc,
                   "helixforge_release": provenance["helixforge_release"], "helixforge_commit": provenance["helixforge_commit"],
                   "caller": groups[0]["caller"], "caller_version": groups[0]["caller_version"],
                   "samples": samples, "annotation_policies": policies, "limitations": entry["interpretation_warning"],
                   "consensus_segments": len(annotated), "intergenic_segments": sum(r["category"] == "intergenic" for r in annotated),
                   "gene_associations": len(associations), "genes_with_associations": len({r["gene_id"] for r in associations}),
                   "stage_segments": dict(Counter(contexts[r["source_id"]][0] for r in annotated)),
                   "coordinate_convention": "0-based_half-open", "scope": "all_atlas_genes",
                   "source_files": [{"path": r["path"], "sha256": r["sha256"]} for r in source["files"]]}
        profiles.append(profile)
    for row in metadata:
        for key,value in {"consensus_method":"replicate_support","annotation_provider":"python_interval_v1","min_libraries":2}.items():row.setdefault(key,value)
    for row in evidence:
        for key,value in {"n_single_library":0,"n_multiple_libraries":row["n_peaks"],"n_ambiguous":0,"n_ambiguity_not_assessed":row["n_peaks"]}.items():row.setdefault(key,value)
    for row in loci:
        for key in SOURCE_FIELDS:row.setdefault(key,{"annotation_provider":"python_interval_v1","association_method":"overlap_priority_first_gene","gene_assignment_ambiguous":"not_assessed","alternative_gene_ids":"[]","occupancy_status":"not_assessed"}.get(key,""))
    stage_order = {"miracidium": 0, "sporocyst": 1, "cercaria": 2, "schistosomulum": 3, "juvenile": 4, "adult": 5}
    metadata.sort(key=lambda row: (stage_order.get(row["stage_or_condition"], 99), row["study_id"], row["mark_or_factor"]))
    for profile in profiles:
        profile["stage_segments"] = dict(sorted(profile["stage_segments"].items(), key=lambda item: stage_order.get(item[0], 99)))
    return {"metadata": metadata, "evidence": evidence, "positions": positions, "loci": loci,
            "segments": segments, "profiles": profiles, "inputs": inputs, "library_peaks": library_peaks, "core_segments":core_segments, "regional_support":regional_support, "idr_segments":idr_segments}


def render_study_page(profile: dict) -> str:
    if profile.get("admission_scope") == "published_descriptive_catalog_only":
        from hf_chipseq_catalog import render_catalog_page
        return render_catalog_page(profile)
    """Small standalone study page: no RNA matrices are loaded to inspect QC."""
    esc = lambda value: html.escape(str(value))
    study = profile["study_id"]
    regional = profile.get("consensus_method") == "nfcore_2_1"
    catalog = "Regional candidates" if regional else "Atomic support segments"
    consensus_fact = "Merged regional candidates · at least 1 calling library; shared core is separate" if regional else "Shared atomic bases from at least 2 deposited IP libraries"
    annotation_tip = "HOMER 4.11 preserves its nearest-promoter assignment and native location class. An intergenic location can still have an associated nearest gene. This is proximity, not evidence of regulation. Ambiguity is not assessed in this source execution; that is not equivalent to an unambiguous assignment." if regional else "Native overlap priority and first-gene assignment; the atlas preserves source classes."
    core_download = '<a href="data/chipseq_shared_core.bed" download>Shared core BED</a><a href="data/chipseq_regional_library_support.tsv" download>Per-region library coverage</a>' if regional else ''
    method_line = '<div>Regional catalog: nf-core/chipseq 2.1.0 · HOMER 4.11 · source execution '+esc(profile.get("regional_execution_commit",""))+'</div>' if regional else ''
    support_cards = ''
    if regional:
        groups = profile["publication"]["groups"]
        single = sum(g["single_library_regions"] for g in groups)
        multiple = sum(g["multiple_library_regions"] for g in groups)
        support_cards = f'<div class="card"><small>1 calling library</small><strong>{single}</strong><span>candidate regions</span></div><div class="card"><small>≥2 calling libraries</small><strong>{multiple}</strong><span>candidate regions; not occupancy significance</span></div>'
    qc_rows = []
    for row in profile["samples"]:
        control = row["is_control"].lower() == "true"
        frip = "—" if control else f'{float(row["frip"]) * 100:.3f}%'
        qc_rows.append('<tr>' + ''.join(f'<td>{esc(value)}</td>' for value in (
            row["sample_id"], row["run_accession"], row["condition"], "Input" if control else "HP1 IP",
            row["final_bam_reads"], "—" if control else row["peak_count"], frip, row["qc_status"],
        )) + '</tr>')
    stage_cards = ''.join(
        f'<div class="card"><small>{esc(stage)}</small><strong>{count}</strong><span>{catalog.lower()}</span></div>'
        for stage, count in profile["stage_segments"].items())
    policies = ''.join(f'<tr><td>{esc(p["condition"])}</td><td>{esc(p["mode"])}</td><td>{esc(p["gene_assignment"])}</td><td>{p["promoter_upstream"]} / {p["promoter_downstream"]} bp</td><td>{esc(p["strand_aware"])}</td></tr>' for p in profile["annotation_policies"])
    policy_headers = '<th>Stage</th><th>Provider</th><th>Assignment</th>' if regional else '<th>Stage</th><th>Mode</th><th>Assignment</th><th>Promoter upstream / downstream</th><th>Strand aware</th>'
    if regional:
        policies = ''.join(f'<tr><td>{esc(p["condition"])}</td><td>HOMER 4.11</td><td>Nearest promoter · ambiguity not assessed</td></tr>' for p in profile["annotation_policies"])
    sources = ''.join(f'<tr><td><a href="{profile["repository"]}/blob/{profile["commit"]}/{r["path"]}">{esc(r["path"])}</a></td><td class="hash">{r["sha256"]}</td></tr>' for r in profile["source_files"])
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{study} · HP1 ChIP-seq · HelixForge atlas</title>
<style>*{{box-sizing:border-box}}body{{margin:0;background:#f3f7f9;color:#152e46;font:15px/1.55 system-ui,Arial,sans-serif}}header{{padding:34px max(24px,calc((100vw - 1180px)/2));background:linear-gradient(125deg,#0b2731,#155c57);color:white}}header small{{letter-spacing:.12em}}h1{{margin:8px 0;font-size:34px}}nav{{background:white;border-bottom:1px solid #d5e1e8;padding:14px max(24px,calc((100vw - 1180px)/2));display:flex;gap:24px;flex-wrap:wrap}}a{{color:#176c71}}main{{max-width:1228px;margin:auto;padding:24px}}h2{{font-size:23px}}h3{{margin:0}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(175px,1fr));gap:14px}}.card,.panel{{background:white;border:1px solid #d5e1e8;border-radius:14px;padding:20px}}.card small,.card span{{display:block;color:#61768a}}.card strong{{display:block;font-size:30px;margin:5px 0}}.panel{{margin-top:20px}}.heading{{display:flex;align-items:center;gap:8px;margin-bottom:15px}}.tip{{position:relative}}.tip summary{{display:grid;place-items:center;list-style:none;width:26px;height:26px;cursor:pointer;background:#edf5f6;border:1px solid #afc8cd;border-radius:50%;font:italic bold 14px Georgia}}.tip summary::-webkit-details-marker{{display:none}}.tip>div{{position:absolute;z-index:5;left:0;top:33px;width:min(350px,75vw);padding:15px;background:white;border:1px solid #cbd8e0;border-radius:10px;box-shadow:0 10px 28px #16314335;font-size:13px}}.scroll{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:13px}}td,th{{padding:10px;border-bottom:1px solid #e1e9ee;text-align:left;vertical-align:top}}th{{background:#edf4f7}}.warning{{border-left:4px solid #c69726;background:#fff9e9;padding:14px 18px;border-radius:8px;margin:20px 0}}.facts{{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:14px}}.facts strong{{display:block}}.hash{{font:11px monospace;overflow-wrap:anywhere;min-width:180px}}.downloads{{display:flex;gap:18px;flex-wrap:wrap}}.disclosure>summary{{cursor:pointer;font-weight:650;padding:10px 0}}footer{{text-align:center;padding:25px;color:#61768a}}@media(max-width:600px){{main{{padding:18px}}h1{{font-size:28px}}.tip>div{{left:auto;right:-20px}}}}@media print{{.tip>div{{position:static;display:block;box-shadow:none}}}}</style></head><body>
<header><small>HELIXFORGE · CHIP-SEQ STUDY</small><h1>{study} · HP1</h1><span>Original per-library processing: HelixForge {profile["helixforge_release"]}</span>{method_line}</header>
<nav><a href="chipseq.html?study={study}">ChIP-seq explorer</a><a href="chipseq.html#studies">← Browse ChIP-seq studies</a><a href="atlas.html">RNA-seq atlas</a><a href="atlas.html#genes">Find a gene</a></nav><main>
<div class="cards"><div class="card"><small>IP libraries</small><strong>{profile["qc"]["ip_libraries"]}</strong></div><div class="card"><small>Input libraries</small><strong>{profile["qc"]["input_libraries"]}</strong></div><div class="card"><small>{catalog}</small><strong>{profile["consensus_segments"]}</strong></div><div class="card"><small>Genes with associations</small><strong>{profile["genes_with_associations"]}</strong></div></div>
<div class="warning"><strong>Descriptive candidates · reviewed with limitations.</strong> Stage and strain differ together; no differential-binding comparison is available.</div>
<div class="panel"><div class="heading"><h2>Experimental design &amp; processing</h2><details class="tip"><summary aria-label="Study interpretation limits">i</summary><div>{esc(profile["limitations"])}</div></details></div>
<div class="facts"><div><strong>Cercaria</strong>SmVEN · 3 deposited HP1 IP libraries + 1 Input</div><div><strong>Sporocyst</strong>SmDFO · 3 deposited HP1 IP libraries + 1 Input</div><div><strong>Reference</strong>{esc(profile["reference_id"])}</div><div><strong>Peak caller</strong>{esc(profile["caller"])} {profile["caller_version"]} · broad</div><div><strong>Catalog</strong>{consensus_fact}</div><div><strong>Extension</strong>200 bp model-free operational assumption</div></div></div>
<div class="panel"><div class="heading"><h2>{catalog}</h2><details class="tip"><summary aria-label="About consensus support segments">i</summary><div>Candidate regions and shared core are distinct products. Library membership does not imply shared support at every base or confirmed occupancy. HOMER nearest-promoter assignments preserve intergenic location labels. Downloadable BED coordinates are 0-based half-open.</div></details></div><div class="cards">{stage_cards}{support_cards}<div class="card"><small>Intergenic</small><strong>{profile["intergenic_segments"]}</strong><span>source location class; nearest-gene links preserved</span></div><div class="card"><small>Gene associations</small><strong>{profile["gene_associations"]}</strong><span>native annotations retained</span></div></div></div>
<div class="panel"><div class="heading"><h2>Quality control</h2><details class="tip"><summary aria-label="About HP1 QC">i</summary><div>All six IP libraries retain REVIEW status. FRiP is the fraction of reads in the called per-library peaks. Low enrichment and extension sensitivity are documented in the accepted study review. Source-library independence is not established.</div></details></div><details class="disclosure"><summary>View all 8 libraries</summary><div class="scroll"><table><thead><tr><th>Library</th><th>Run</th><th>Stage</th><th>Type</th><th>Final BAM reads</th><th>Called peaks</th><th>FRiP</th><th>QC</th></tr></thead><tbody>{''.join(qc_rows)}</tbody></table></div></details></div>
<div class="panel"><div class="heading"><h2>Annotation policy</h2><details class="tip"><summary aria-label="About native gene assignment">i</summary><div>{annotation_tip}</div></details></div><details class="disclosure"><summary>View native annotation settings</summary><div class="scroll"><table><thead><tr>{policy_headers}</tr></thead><tbody>{policies}</tbody></table></div></details></div>
<div class="panel"><h2>Provenance &amp; downloads</h2><div class="downloads"><a href="data/chipseq_{study}.json" download>Study design, QC &amp; provenance</a><a href="data/chipseq_{study}_qc.tsv" download>Per-library QC</a><a href="data/chipseq_regional_candidates.tsv" download>All regional candidates</a>{core_download}<a href="data/chipseq_gene_peak_loci.tsv" download>All nearest-gene associations</a><a href="{profile["repository"]}/tree/{profile["commit"]}/results/chipseq/{study}">Published HF package</a></div><p>Source commit: <code>{profile["commit"]}</code></p><details class="disclosure"><summary>Frozen source files and SHA-256 checksums</summary><div class="scroll"><table><thead><tr><th>Source file</th><th>SHA-256</th></tr></thead><tbody>{sources}</tbody></table></div></details></div>
</main><footer>HelixForge-SMansoni · study-specific candidate evidence</footer></body></html>'''


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", required=True, type=Path)
    parser.add_argument("--revision", required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    freeze_source(args.repository, args.revision, args.study, args.output)
