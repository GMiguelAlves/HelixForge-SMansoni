"""Read published nf-core/HOMER catalogs without recomputing associations."""
import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict

SOURCE_FIELDS = ["annotation_provider", "association_method", "distance_to_tss", "calling_libraries",
                 "total_library_count", "library_support_category", "support_core_bp",
                 "gene_assignment_ambiguous", "gene_assignment_ambiguity_reason",
                 "gene_assignment_ambiguity_scope", "alternative_gene_ids", "occupancy_status"]


def rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bed(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip().split("\t") for line in handle if line.strip()]


def load_regional_study(root, result, entry, source, config, annotations):
    publication = json.loads((result / "manifests/regional_publication.json").read_text(encoding="utf-8"))
    qc = json.loads((result / "qc/qc_summary.json").read_text(encoding="utf-8"))
    if (qc["atlas_inclusion"], qc["technical_execution"], qc["processing_status"]) != ("ACCEPTED", "PASS", "COMPLETE"):
        raise ValueError("Regional study has not been accepted")
    if any(publication[k] != value for k, value in {"status": "PUBLISHED", "publication_checks": "PASS",
            "source_workflow_validation": "PASS", "source_direct_equivalence": "PASS",
            "published_representation": "nfcore_2_1_min1_homer_4_11"}.items()):
        raise ValueError("Regional publication failed its declared admission gate")
    reference = json.loads((root / "config/reference/chipseq_chromosomes_wbps19.json").read_text(encoding="utf-8"))
    chroms = reference["chromosomes"]
    bundle = json.loads((result / "manifests/reference_bundle.manifest.json").read_text(encoding="utf-8"))
    if bundle["genome_id"] != config["reference_id"]:raise ValueError("Incompatible regional reference bundle")
    caller_versions = set()
    def interval(chrom, start, end):
        if chrom not in chroms or not 0 <= int(start) < int(end) <= chroms[chrom]:
            raise ValueError("Regional interval/reference mismatch")
    samples = rows(result / "qc/qc_sample_summary.tsv")
    sample_map = {r["run_accession"]: r for r in samples}
    output = {key: [] for key in ("metadata", "evidence", "positions", "loci", "segments", "profiles", "library_peaks", "core_segments", "regional_support")}
    summaries, positions = {}, defaultdict(set)
    study = entry["id"]
    policies, contexts = [], []
    category_counts, stage_counts = Counter(), Counter()
    for group in publication["groups"]:
        stage, factor = group["condition"], group["target"]
        context = f"{study}|{stage}|{factor}"
        directory = result / "broad/chipseq/consensus" / stage
        manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
        if manifest["reference_id"] != config["reference_id"] or manifest["reference_sha256"] != reference["reference_fasta_sha256"]:
            raise ValueError("Incompatible regional reference")
        if manifest["provider"] != "nfcore_chipseq_2.1.0" or manifest["annotation_provider"] != "homer_4_11" or manifest["regional_min_libraries"] != 1 or manifest["core_min_libraries"] != 2:
            raise ValueError("Unexpected regional method or support policy")
        for artifact in manifest["artifacts"].values():
            path = result / artifact["path"]
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != artifact["sha256"]:
                raise ValueError("Regional public artifact checksum mismatch")
        if manifest["annotation_sha256"] != bundle["artifacts"]["annotation"]["sha256"]:raise ValueError("Regional GTF checksum differs from reference bundle")
        declared = manifest["library_metadata"]
        libs = {r["record_id"]: r for r in declared}
        if len(libs) != group["total_library_count"] or len(libs) != len(declared):
            raise ValueError("Invalid regional library identities")
        for rid, library in libs.items():
            sample = sample_map[rid]
            if sample["sample_id"] != library["sample_id"] or sample["condition"] != stage or sample["is_control"] != "false":
                raise ValueError("Regional library metadata differs from QC")
            peak_path = result / f"broad/080-peak-calling/{rid}.{factor}.broad.macs3.peak_calling/peaks.broadPeak"
            caller_manifest = json.loads((peak_path.parent / "manifest.json").read_text(encoding="utf-8"))
            if caller_manifest["parameters"]["reference_sha256"] != reference["reference_fasta_sha256"] or caller_manifest["control_record_id"] != library["control_record_id"] or hashlib.sha256(peak_path.read_bytes()).hexdigest() != caller_manifest["artifacts"]["peaks"]["sha256"]:raise ValueError("Native caller reference/control/checksum mismatch")
            caller_versions.add((caller_manifest["caller"], caller_manifest["caller_version"]))
            calls = []
            for line in peak_path.read_text(encoding="utf-8").splitlines():
                if not line.strip():
                    continue
                p = line.split("\t")
                if len(p) != 9:
                    raise ValueError("Invalid native broadPeak columns")
                interval(p[0], p[1], p[2])
                calls.append({"study_id": study, "stage_or_condition": stage, "mark_or_factor": factor,
                    "sample_id": sample["sample_id"], "run_accession": rid, "peak_id": p[3], "chrom": p[0],
                    "start": int(p[1]), "end": int(p[2]), "signal_value": float(p[6]), "p_value": float(p[7]), "q_value": float(p[8])})
            if len(calls) != int(sample["peak_count"]):
                raise ValueError("Native regional-library calls differ from QC")
            output["library_peaks"].extend(calls)
        native_rows = rows(directory / "consolidated_peaks.tsv.gz")
        native = {r["peak_id"]: r for r in native_rows}
        merged = []
        for call in sorted((r for r in output["library_peaks"] if r["stage_or_condition"]==stage and r["mark_or_factor"]==factor),key=lambda r:(r["chrom"],r["start"],r["end"])):
            if merged and call["chrom"]==merged[-1][0] and call["start"]<=merged[-1][2]:
                merged[-1][2]=max(merged[-1][2],call["end"])
            else:merged.append([call["chrom"],call["start"],call["end"]])
        if {tuple(r) for r in merged}!={(r["chrom"],int(r["start"]),int(r["end"])) for r in native_rows}:
            raise ValueError("Minimum-one regional candidates do not preserve the native called-interval union")

        annotation = rows(result / f"broad/chipseq/peak_annotation/{stage}/annotated_peaks.tsv.gz")
        assoc = rows(result / f"broad/chipseq/peak_annotation/{stage}/peak_gene_associations.tsv.gz")
        annotated = {r["peak_id"]: r for r in annotation}
        if len(native) != len(native_rows) or len(annotated) != len(annotation) or native.keys() != annotated.keys():
            raise ValueError("Regional annotation does not conserve unique candidate IDs")
        support_by_peak = defaultdict(dict)
        for row in rows(directory / "support_by_library.tsv.gz"):
            row = {k:v for k,v in row.items() if k != "library_support_json"}
            pid, rid = row["peak_id"], row["record_id"]
            if pid not in native or rid not in libs or rid in support_by_peak[pid] or row["sample_id"] != libs[rid]["sample_id"]:
                raise ValueError("Duplicate/incompatible region-library evidence")
            length = int(native[pid]["end"]) - int(native[pid]["start"])
            covered, core, calls = int(row["covered_bp"]), int(row["region_core_bp"]), int(row["called_peak_count"])
            fraction = float(row["covered_fraction"])
            if int(row["region_bp"]) != length or not 0 <= covered <= length or not 0 <= core <= length or calls < 0 or (calls > 0) != (covered > 0) or not math.isfinite(fraction) or abs(fraction - covered / length) > 1e-6 or row["occupancy_status"] != "not_assessed":
                raise ValueError("Invalid regional coverage evidence")
            support_by_peak[pid][rid] = row
            output["regional_support"].append({"study_id": study, "context_id": context, **row})
        if sum(map(len, support_by_peak.values())) != len(native) * len(libs):
            raise ValueError("Regional evidence omits libraries with zero support")
        for rid in libs:
            if sum(int(details[rid]["called_peak_count"]) for details in support_by_peak.values()) != int(sample_map[rid]["peak_count"]):
                raise ValueError("Regional support does not conserve every original library call")
        source_associations = set()
        for pid, row in annotated.items():
            region = native[pid]
            interval(row["chrom"], row["start"], row["end"])
            if any(row[k] != region[k] for k in ("chrom", "start", "end")):
                raise ValueError("Regional annotation changed native coordinates")
            details = support_by_peak[pid]
            calling = sorted(rid for rid, d in details.items() if int(d["called_peak_count"]) > 0)
            category = "single_library" if len(calling) == 1 else "all_libraries" if len(calling) == len(libs) else "multiple_libraries"
            classification = "candidate" if len(calling) == 1 else "candidate_with_library_support"
            expected = {"calling_libraries": ",".join(calling), "calling_samples": ",".join(libs[rid]["sample_id"] for rid in calling),
                        "calling_library_count": str(len(calling)), "total_library_count": str(len(libs)), "library_support_category": category}
            if not calling or any(region[k] != v or row[k] != v for k, v in expected.items()) or int(region["support"]) != len(calling) or region["classification"] != classification or row["evidence_classification"] != classification:
                raise ValueError("Regional library summary/classification differs from evidence")
            decoded = json.loads(row["library_support_json"])
            if {d["record_id"]: d for d in decoded} != details or len(decoded) != len(libs):
                raise ValueError("Annotation lost or changed per-library coverage evidence")
            cores = {int(d["region_core_bp"]) for d in details.values()}
            if len(cores) != 1 or int(row["support_core_bp"]) not in cores:
                raise ValueError("Shared core differs across evidence rows")
            if row["annotation_provider"] != "homer_4_11" or row["association_method"] != "homer_nearest_promoter" or row["occupancy_status"] != "not_assessed":
                raise ValueError("Regional annotation method/status changed")
            flag = row["gene_assignment_ambiguous"]
            alternatives = json.loads(row["alternative_gene_ids"])
            if manifest["gene_assignment_ambiguity"] == "not_assessed" and flag != "not_assessed":raise ValueError("Regional ambiguity state differs from source assessment")
            if flag not in {"true", "false", "not_assessed"} or any(g not in annotations for g in alternatives) or flag == "true" and not alternatives:
                raise ValueError("Invalid gene-assignment ambiguity evidence")
            selected = list(filter(None, row["gene_ids"].split(';')))
            if row["gene_count"] != str(len(selected)) or any(g not in annotations for g in selected):
                raise ValueError("Regional assigned gene is outside the reference")
            source_associations.update((pid, g) for g in selected)
            extras = {k: row[k] for k in SOURCE_FIELDS}
            extras.update(total_library_count=int(row["total_library_count"]), support_core_bp=int(row["support_core_bp"]))
            output["segments"].append({"study_id": study, "context_id": context, "stage_or_condition": stage, "mark_or_factor": factor,
                **{k: row[k] for k in ("peak_id", "chrom", "start", "end", "category", "gene_ids")},
                "support": len(calling), "support_libraries": row["calling_libraries"], "consensus_method": "nfcore_2_1", **extras})
            category_counts[row["category"]] += 1
            stage_counts[stage] += 1
        seen = set()
        for row in assoc:
            pid, gene = row["peak_id"], row["gene_id"]
            if (pid, gene) in seen or (pid, gene) not in source_associations:
                raise ValueError("Duplicate/inconsistent native regional association")
            seen.add((pid, gene))
            peak = annotated[pid]
            if annotations[gene]["chromosome"] != peak["chrom"] or any(row[k] != peak[k] for k in SOURCE_FIELDS + ["category", "library_support_json"]):
                raise ValueError("Native region/association evidence differs")
            category = row["category"]
            position = "promoter_TSS" if category == "promoter" else "gene_body" if category in {"exon", "intron", "gene"} else "distal_intergenic" if category == "intergenic" else "other_annotated"
            summary = summaries.setdefault((gene, context), {"study_id": study, "context_id": context, "gene_id": gene,
                "functional_name": annotations[gene]["functional_name"], "mark_or_factor": factor, "stage_or_condition": stage,
                "n_peaks": 0, "n_promoter_peaks": 0, "peak_locations": set(), "regulatory_class": "candidate_nearest_promoter_association",
                "expected_effect": "not_inferred_from_occupancy", "n_single_library": 0, "n_multiple_libraries": 0,
                "n_ambiguous": 0, "n_ambiguity_not_assessed": 0})
            summary["n_peaks"] += 1
            summary["n_promoter_peaks"] += int(category == "promoter")
            summary["peak_locations"].add(category)
            summary["n_single_library"] += int(peak["calling_library_count"]) == 1
            summary["n_multiple_libraries"] += int(peak["calling_library_count"]) >= 2
            summary["n_ambiguous"] += peak["gene_assignment_ambiguous"] == "true"
            summary["n_ambiguity_not_assessed"] += peak["gene_assignment_ambiguous"] == "not_assessed"
            positions[(gene, context, stage, factor, position)].add(pid)
            output["loci"].append({"study_id": study, "context_id": context, "gene_id": gene, "stage_or_condition": stage,
                "mark_or_factor": factor, "peak_id": pid, "peak_chrom": peak["chrom"], "peak_start": int(peak["start"]),
                "peak_end": int(peak["end"]), "peak_midpoint": (int(peak["start"]) + int(peak["end"])) // 2,
                "peak_location": category, "promoter_flag": str(category == "promoter").lower(), "position_class": position,
                "support": int(native[pid]["support"]), **{k: peak[k] for k in SOURCE_FIELDS}})
        if seen != source_associations or len(native) != group["regions"] or len({r["gene_id"] for r in assoc}) != group["associated_gene_ids"] or sum(int(r["end"])-int(r["start"]) for r in native.values()) != group["region_bp"]:
            raise ValueError("Published regional counts/coverage/associations not conserved")
        core_rows = bed(directory / "support_core.bed.gz")
        for chrom, start, end, *_ in core_rows:
            interval(chrom, start, end)
            if not any(r["chrom"] == chrom and int(r["start"]) <= int(start) < int(end) <= int(r["end"]) for r in native.values()):
                raise ValueError("Shared core extends outside regional candidates")
            output["core_segments"].append({"study_id": study, "stage_or_condition": stage, "mark_or_factor": factor, "chrom": chrom, "start": int(start), "end": int(end)})
        if sum(int(r[2])-int(r[1]) for r in core_rows) != group["core_bp"]:
            raise ValueError("Shared-core coverage differs from publication")
        output["metadata"].append({"study_id": study, "context_id": context, "mark_or_factor": factor, "stage_or_condition": stage,
            "n_samples": len(libs), "n_controls": len({r["control_record_id"] for r in declared}), "replicates": ";".join(r["sample_id"] for r in declared),
            "batches": "not_established", "qc_status": qc["qc_classification"], "consensus_method": "nfcore_2_1", "annotation_provider": "homer_4_11", "min_libraries": 1})
        contexts.append({"study": study, "stage": stage, "factor": factor, "method": "nfcore_2_1", "min_libraries": 1,
            "libraries": [{"id": rid, "sample": libs[rid]["sample_id"], "label": f"Library {i+1}"} for i, rid in enumerate(sorted(libs))]})
        policies.append({"condition": stage, "provider": "homer_4_11", "mode": "homer_native", "gene_assignment": "homer_nearest_promoter", "strand_aware": "provider_native", "promoter_upstream": "provider_native", "promoter_downstream": "provider_native"})
    output["evidence"] = [{**r, "peak_locations": ";".join(sorted(r["peak_locations"]))} for r in summaries.values()]
    output["positions"] = [{"study_id": study, "context_id": context, "gene_id": gene, "stage_or_condition": stage, "mark_or_factor": factor,
        "position_class": position, "n_unique_peaks": len(ids)} for (gene, context, stage, factor, position), ids in positions.items()]
    provenance = json.loads((result.parents[2] / f"provenance/chipseq/{study}/execution_provenance.json").read_text(encoding="utf-8"))
    if len(caller_versions) != 1:raise ValueError("Mixed peak caller versions in regional catalog")
    caller, caller_version = next(iter(caller_versions))
    profile = {"study_id": study, "reference_id": config["reference_id"], "target": publication["groups"][0]["target"],
        "repository": source["repository"], "commit": source["commit"], "qc": qc, "helixforge_release": provenance["helixforge_release"],
        "helixforge_commit": provenance["helixforge_commit"], "regional_execution_commit": manifest["helixforge_execution_commit"],
        "caller": caller, "caller_version": caller_version, "consensus_method": "nfcore_2_1", "annotation_provider": "homer_4_11",
        "regional_min_libraries": 1, "core_min_libraries": 2, "contexts": contexts, "samples": samples, "annotation_policies": policies,
        "limitations": entry["interpretation_warning"], "consensus_segments": len(output["segments"]), "intergenic_segments": category_counts["intergenic"],
        "gene_associations": len(output["loci"]), "genes_with_associations": len({r["gene_id"] for r in output["loci"]}), "stage_segments": dict(stage_counts),
        "coordinate_convention": "0-based_half-open", "scope": "all_atlas_genes", "publication": publication,
        "source_files": [{"path": r["path"], "sha256": r["sha256"]} for r in source["files"]]}
    output["profiles"].append(profile)
    return output
