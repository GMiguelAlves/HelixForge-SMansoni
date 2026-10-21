#!/usr/bin/env python3
"""Validate the committed RNA-seq atlas bundle without optional dependencies."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


EXPECTED = {
    "studies": ["PRJNA602528", "PRJNA597909", "PRJEB14695", "PRJEB32839", "E-MTAB-451", "PRJEB3190", "E-ERAD-478", "PRJNA704646", "PRJNA1150591", "PRJNA312093", "PRJNA1177852", "PRJNA428470"],
    "biological_samples": 268,
    "technical_runs": 535,
    "genes": 9914,
    "contrasts": 71,
}
EXPECTED_EXTERNAL_TYPES = {
    "PUBLISHED_STAGE_ENRICHMENT": 1797,
    "PUBLISHED_GENE_MODEL_EXTENSION": 6819,
    "PUBLISHED_MODEL_RESTRUCTURING": 668,
    "PUTATIVE_NOVEL_GENE_CLAIM": 232,
}
EXPECTED_LNCRNA_CLASSES = {"LINC_RNA": 7029, "ANTISENSE_LNCRNA": 402}
EXPECTED_ALMEIDA_TYPES = {
    "PUBLISHED_ADULT_MALE_TRANSCRIPTION_SUPPORT": 9928,
    "PUBLISHED_FIRST_ADULT_MALE_TRANSCRIPTION_SUPPORT": 990,
    "PUBLISHED_NOVEL_TRANSCRIBED_CONTIG": 989,
}
FORBIDDEN = (
    "/ho" + "me/",
    "/scr" + "atch/",
    "srv-" + "slurm",
    "@bio." + "ib.unicamp.br",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() and path.with_suffix(path.suffix+".gz").is_file():path=path.with_suffix(path.suffix+".gz")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(root: Path, atlas: Path | None = None) -> list[str]:
    errors: list[str] = []
    atlas = atlas or root / "results/atlas/rnaseq"
    manifest_path = atlas / "manifest.json"
    if not manifest_path.is_file():
        return ["atlas manifest is missing"]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for field, expected in EXPECTED.items():
        if manifest.get(field) != expected:
            errors.append(f"manifest {field}: expected {expected!r}, observed {manifest.get(field)!r}")
    if manifest.get("status") != "PASS":
        errors.append("atlas status is not PASS")
    display = manifest.get("locus_display", {})
    expected_display = json.loads((root / "analyses/atlas/atlas_config.json").read_text(encoding="utf-8"))["locus_display"]
    if display != expected_display or display.get("reclassify_source_peaks") is not False or display.get("basis") != "STRAND_AWARE_DISPLAY_WINDOWS_NATIVE_SOURCE_CLASS_PRESERVED":
        errors.append("locus reference-window policy is missing or changes source peak classification")
    if manifest.get("batch_effect_assessment") != "DEFERRED":
        errors.append("batch assessment must remain DEFERRED in atlas v1")
    if manifest.get("global_pca_interpretation") != "EXPLORATORY_NOT_A_BATCH_TEST":
        errors.append("global PCA interpretation boundary is missing")
    if manifest.get("chipseq_integration_status") != "HF_PROCESSED_REVIEWED":
        errors.append("reviewed native HF ChIP-seq status is missing")
    if manifest.get("chipseq_scope") != "all_atlas_genes":
        errors.append("ChIP evidence is not available across the full atlas gene universe")

    for record in manifest.get("inputs", []):
        path = root / record["path"]
        if not path.is_file():
            errors.append(f"missing input: {record['path']}")
            continue
        if sha256(path) != record["sha256"]:
            errors.append(f"input checksum mismatch: {record['path']}")

    for record in manifest.get("outputs", []):
        path = atlas / record["path"]
        if not path.is_file():
            errors.append(f"missing output: {record['path']}")
            continue
        if path.stat().st_size != record["bytes"]:
            errors.append(f"size mismatch: {record['path']}")
        if sha256(path) != record["sha256"]:
            errors.append(f"checksum mismatch: {record['path']}")

    samples = read_tsv(atlas / "data/sample_inventory.tsv")
    if len(samples) != EXPECTED["biological_samples"]:
        errors.append("sample inventory cardinality mismatch")
    if sum(int(row["technical_runs"]) for row in samples) != EXPECTED["technical_runs"]:
        errors.append("technical-run total mismatch")
    if {row["reference_id"] for row in samples} != {"Schistosoma_mansoni_SM_V10_WBPS19"}:
        errors.append("sample inventory mixes reference identities")
    if any(not row["batch"] for row in samples):
        errors.append("atlas contains empty batch metadata")
    if {row["study"] for row in samples} != set(EXPECTED["studies"]):
        errors.append("external literature was incorrectly added as a bulk study")

    contrasts = read_tsv(atlas / "data/contrast_catalog.tsv")
    if len(contrasts) != EXPECTED["contrasts"]:
        errors.append("contrast catalog cardinality mismatch")
    if {row["study"] for row in contrasts} != {"PRJNA597909", "PRJEB14695", "PRJEB32839", "E-MTAB-451", "PRJEB3190", "E-ERAD-478", "PRJNA1150591", "PRJNA312093", "PRJNA428470"}:
        errors.append("unexpected studies in differential-expression catalog")

    allowed = {row["__rowname"] for row in read_tsv(root / "results/PRJNA312093/restricted_refit/included_samples.tsv")}
    observed = {row["sample_key"] for row in samples if row["study"] == "PRJNA312093"}
    if len(allowed) != 12 or observed != allowed:
        errors.append("PRJNA312093 expression does not match the 12 approved libraries")
    for study, n_samples, n_contrasts in (("PRJNA312093", 12, 3), ("PRJNA1177852", 8, 0), ("PRJNA428470", 6, 1)):
        rows = [row for row in contrasts if row["study"] == study]
        if sum(row["study"] == study for row in samples) != n_samples or len(rows) != n_contrasts:
            errors.append(f"admission scope differs for {study}")
        if any(not row.get("interpretation_warning") for row in rows):
            errors.append(f"contrast warning missing for {study}")
    refit = json.loads((root / "results/PRJNA312093/restricted_refit/manifest.json").read_text(encoding="utf-8"))
    refit_rows = [row for row in contrasts if row["study"] == "PRJNA312093"]
    if [int(row["significant_genes"]) for row in refit_rows] != [290, 1045, 437]:
        errors.append("PRJNA312093 significant-gene counts differ from the reviewed refit")
    if any("adult" in row["contrast_key"] or "cercaria" in row["contrast_key"] for row in refit_rows):
        errors.append("non-evaluable PRJNA312093 contrasts were included")
    # Check every exported refit effect and adjusted p-value against its frozen source.
    for filename, source_field in (("gene_contrast_log2fc.tsv", "log2FoldChange"), ("gene_contrast_padj.tsv", "padj")):
        exported = {row["gene_id"]: row for row in read_tsv(atlas / "data" / filename)}
        for record in refit["contrasts"]:
            source = root / "results/PRJNA312093/restricted_refit" / record["path"]
            if sha256(source) != record["sha256"]:
                errors.append(f"frozen refit checksum differs: {record['path']}")
                continue
            key = "PRJNA312093:" + source.stem
            for row in read_tsv(source):
                value = row[source_field]
                expected = "NA" if value in ("NA", "NaN", "") else f"{float(value):.8g}"
                if exported[row["gene_id"]][key] != expected:
                    errors.append(f"refit statistic differs: {filename}, {row['gene_id']}, {key}")
                    break


    annotations = read_tsv(atlas / "data/gene_annotations.tsv")
    annotation_ids = {row["gene_id"] for row in annotations}
    reference_ids = {row["gene_id"] for row in read_tsv(root / "config/reference/gene_annotations_wbps19.tsv")}
    model_source = json.loads((root / "config/reference/transcript_models_wbps19.json").read_text(encoding="utf-8"))
    model_export = json.loads((atlas / "data/gene_models.json").read_text(encoding="utf-8"))
    expected_models = {gene: model_source["genes"].get(gene, []) for gene in annotation_ids}
    if model_export.get("genes") != expected_models or model_export.get("coordinate_system") != "1-based-inclusive":
        errors.append("exported exon structures differ from the frozen annotation")
    if model_export.get("source_sha256") != manifest.get("transcript_model_source_sha256"):
        errors.append("transcript structure source checksum differs")
    if manifest.get("genes_with_transcript_models") != sum(bool(models) for models in expected_models.values()):
        errors.append("transcript structure coverage count differs")
    if manifest.get("transcript_models") != sum(map(len, expected_models.values())):
        errors.append("transcript count differs from the frozen annotation")
    for annotation in annotations:
        if not expected_models[annotation["gene_id"]] and annotation["biotype"] != "pseudogene":
            errors.append(f"exon structure missing for non-pseudogene: {annotation['gene_id']}")
    if len(annotations) != EXPECTED["genes"]:
        errors.append("gene annotation cardinality mismatch")
    if any(not row["functional_name"] for row in annotations):
        errors.append("gene annotations contain empty functional names")
    external = read_tsv(atlas / "data/external_annotation_evidence.tsv")
    markers = read_tsv(atlas / "data/cell_type_marker_evidence.tsv")
    if not external or not markers:
        errors.append("published external evidence is missing")
    if {row["source_study"] for row in external} != {"PRJNA294789"}:
        errors.append("Anderson evidence has an unexpected source identity")
    if {row["source_study"] for row in markers} != {"PRJNA395457"}:
        errors.append("Wang evidence has an unexpected source identity")
    external_types = {kind: sum(row["evidence_type"] == kind for row in external)
                      for kind in {row["evidence_type"] for row in external}}
    if external_types != EXPECTED_EXTERNAL_TYPES:
        errors.append(f"published Anderson source counts differ: {external_types}")
    if len(markers) != 154 or any(not all(f"{cell} TPM=" in row["reported_statistic"] for cell in ("φ", "δ", "κ")) for row in markers):
        errors.append("published Wang class-dependent transcript table is incomplete")
    source_hashes = {entry["sha256"] for entry in manifest.get("inputs", [])
                     if "local_inputs/external_rna/" in entry["path"]}
    for table in (external, markers):
        for row in table:
            if not row["source_publication"] or not row["source_location"] or not row["source_record_id"]:
                errors.append("published evidence lacks source provenance")
                break
            if row["source_checksum"] not in source_hashes:
                errors.append("published evidence checksum is not recorded among source inputs")
                break
            if row["mapping_status"] == "MAPPED":
                if row["current_gene_id"] not in annotation_ids:
                    errors.append("mapped external evidence is outside the atlas gene universe")
                    break
            elif row["mapping_status"] in {"AMBIGUOUS", "UNMAPPED"}:
                if row["current_gene_id"]:
                    errors.append("unresolved external evidence was assigned a gene")
                    break
            else:
                errors.append("unknown external evidence mapping status")
                break
    if manifest.get("external_annotation_records") != len(external) or manifest.get("cell_marker_records") != len(markers):
        errors.append("external evidence manifest cardinalities do not match the tables")
    if manifest.get("external_annotation_mapped") != sum(row["mapping_status"] == "MAPPED" for row in external):
        errors.append("Anderson mapping count differs from the manifest")
    if manifest.get("cell_marker_mapped") != sum(row["mapping_status"] == "MAPPED" for row in markers):
        errors.append("Wang mapping count differs from the manifest")

    lncrnas = read_tsv(atlas / "data/lncrna_annotation_evidence.tsv")
    legacy = read_tsv(atlas / "data/legacy_transcript_evidence.tsv")
    lncrna_classes = {kind: sum(row["lncrna_class"] == kind for row in lncrnas)
                      for kind in {row["lncrna_class"] for row in lncrnas}}
    if lncrna_classes != EXPECTED_LNCRNA_CLASSES:
        errors.append(f"Vasconcelos lncRNA source counts differ: {lncrna_classes}")
    if sum(row["robust_candidate"] == "YES" for row in lncrnas) != 181:
        errors.append("Vasconcelos robust-candidate count differs")
    if sum(row["stage_differential_expression"] == "YES" for row in lncrnas) != 916:
        errors.append("Vasconcelos stage-DE count differs")
    if any(not row["source_chromosome"] or not row["source_start"] or not row["source_end"] for row in lncrnas):
        errors.append("Vasconcelos catalogue contains an incomplete source coordinate")
    for row in lncrnas:
        mapped = [gene for gene in row["associated_current_gene_id"].split(";") if gene]
        if any(gene not in annotation_ids for gene in mapped):
            errors.append("Vasconcelos antisense association is outside the atlas gene universe")
            break
        if row["association_mapping_status"] in {"UNMAPPED", "AMBIGUOUS", "NOT_APPLICABLE"} and mapped:
            errors.append("unresolved Vasconcelos association was assigned a current gene")
            break
    almeida_types = {kind: sum(row["evidence_type"] == kind for row in legacy)
                     for kind in {row["evidence_type"] for row in legacy}}
    if almeida_types != EXPECTED_ALMEIDA_TYPES:
        errors.append(f"Almeida source counts differ: {almeida_types}")
    if any(row["current_gene_id"] and row["current_gene_id"] not in annotation_ids for row in legacy):
        errors.append("Almeida evidence is linked outside the atlas gene universe")
    if any(row["evidence_type"] == "PUBLISHED_NOVEL_TRANSCRIBED_CONTIG" and row["current_gene_id"] for row in legacy):
        errors.append("an Almeida historical contig was assigned to a current gene")
    literature_hashes = {entry["sha256"] for entry in manifest.get("inputs", [])
                         if "local_inputs/literature_rna/" in entry["path"]}
    if len(literature_hashes) != 6:
        errors.append("the six Vasconcelos/Almeida source tables are not all checksum-recorded")
    for row in lncrnas:
        if row["catalog_source_checksum"] not in literature_hashes:
            errors.append("Vasconcelos catalogue checksum is absent from manifest inputs")
            break
        if row["stage_de_source_checksum"] and row["stage_de_source_checksum"] not in literature_hashes:
            errors.append("Vasconcelos stage-DE checksum is absent from manifest inputs")
            break
        if row["robust_source_checksum"] and row["robust_source_checksum"] not in literature_hashes:
            errors.append("Vasconcelos robust-candidate checksum is absent from manifest inputs")
            break
    if any(row["source_checksum"] not in literature_hashes for row in legacy):
        errors.append("Almeida evidence checksum is absent from manifest inputs")
    expected_lncrna_links = sum(len([gene for gene in row["associated_current_gene_id"].split(";") if gene]) for row in lncrnas)
    expected_almeida_mapped = sum(bool(row["current_gene_id"]) for row in legacy)
    if manifest.get("vasconcelos_lncrna_records") != len(lncrnas):
        errors.append("Vasconcelos manifest cardinality differs from the evidence table")
    if manifest.get("vasconcelos_current_gene_links") != expected_lncrna_links:
        errors.append("Vasconcelos manifest gene-link count differs")
    if manifest.get("almeida_evidence_records") != len(legacy) or manifest.get("almeida_mapped_records") != expected_almeida_mapped:
        errors.append("Almeida manifest cardinality differs from the evidence table")
    selected = read_tsv(atlas / "data/selected_genes.tsv")
    config = json.loads((root / "analyses/atlas/atlas_config.json").read_text(encoding="utf-8"))
    collection_path = atlas / "data/gene_collections.json"
    if not collection_path.is_file():
        errors.append("shared gene collection catalog is missing")
    else:
        collections = json.loads(collection_path.read_text(encoding="utf-8"))
        reference_genes = {row["gene_id"] for row in read_tsv(root / config["gene_annotations"])}
        catalog_genes = {row[0] for row in collections["genes"]}
        machinery_genes = {row["resolved_gene_id"] for row in selected if row["resolved_gene_id"]}
        if collections["reference_id"] != manifest["reference_id"] or catalog_genes != reference_genes:
            errors.append("collection reference catalog differs from the frozen annotation")
        if len(collections["genes"]) != len(catalog_genes) or sum(row[1] for row in collections["genes"]) != EXPECTED["genes"]:
            errors.append("collection catalog duplicates genes or changes RNA availability")
        if set(collections["presets"]["machinery"]["ids"]) != machinery_genes:
            errors.append("machinery collection differs from admitted resolved candidates")
        if {gene for ids in collections["groups"].values() for gene in ids} != machinery_genes:
            errors.append("functional collection shortcuts differ from the machinery")
        if any(gene not in reference_genes for ids in collections["aliases"].values() for gene in ids):
            errors.append("collection alias targets a gene outside the reference")
    candidate_source = root / config["selected_genes"]
    curated = read_tsv(candidate_source)
    exported_source = atlas / "data/histone_modifier_gene_set_v2.tsv"
    if not exported_source.is_file() or sha256(exported_source) != sha256(candidate_source):
        errors.append("full candidate curation export differs from its source")
    if len({row["record_id"] for row in curated}) != len(curated):
        errors.append("curated candidate record IDs are not unique")
    excluded_gene_ids = {row["input_gene_id"] for row in curated if row["scope_tier"] == "EXCLUDED_REVIEW"}
    active = [row for row in curated if row["input_gene_id"] not in excluded_gene_ids]
    by_input: dict[str, list[dict[str, str]]] = {}
    for row in active:
        by_input.setdefault(row["input_gene_id"], []).append(row)
    if {row["gene_id"] for row in selected} != set(by_input) or len(selected) != len(by_input):
        errors.append("selected genes differ from active curated candidates")
    for row in selected:
        source_rows = by_input.get(row["gene_id"], [])
        if not source_rows:
            continue
        if {item["sm_v10_gene_id"] for item in source_rows} != {row["resolved_gene_id"]}:
            errors.append(f"curated mapping differs for {row['gene_id']}")
        if {item["mapping_status"] for item in source_rows} != {row["resolution_status"]}:
            errors.append(f"curated mapping status differs for {row['gene_id']}")
        if {item["record_id"] for item in source_rows} != set(row["record_ids"].split(";")):
            errors.append(f"curated evidence records differ for {row['gene_id']}")
        if row["resolved_gene_id"] not in annotation_ids:
            errors.append(f"selected gene absent from atlas reference: {row['resolved_gene_id']}")
    resolution_counts: dict[str, int] = {}
    for row in selected:
        resolution_counts[row["resolution_status"]] = resolution_counts.get(row["resolution_status"], 0) + 1
    expected_resolution: dict[str, int] = {}
    for source_rows in by_input.values():
        status = source_rows[0]["mapping_status"]
        expected_resolution[status] = expected_resolution.get(status, 0) + 1
    if resolution_counts != expected_resolution:
        errors.append(
            f"candidate ID resolution differs: expected {expected_resolution}, observed {resolution_counts}"
        )
    if manifest.get("selected_source_records") != len(curated) or manifest.get("selected_excluded_review_records") != sum(row["scope_tier"] == "EXCLUDED_REVIEW" for row in curated):
        errors.append("candidate source/exclusion manifest counts differ")
    if manifest.get("selected_mapping_statuses") != expected_resolution:
        errors.append("candidate mapping manifest differs from source")
    if manifest.get("selected_unique_genes") != len(by_input) or manifest.get("selected_unique_genes_present") != len({row["resolved_gene_id"] for row in selected}):
        errors.append("candidate unique-gene manifest counts differ")

    for filename in ("gene_contrast_log2fc.tsv", "gene_contrast_padj.tsv"):
        path = atlas / "data" / filename
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            row_count = sum(1 for _ in reader)
        if len(header) - 1 != EXPECTED["contrasts"] or row_count != EXPECTED["genes"]:
            errors.append(f"unexpected dimensions: {filename}")

    chip_inventory = read_tsv(atlas / "data/chipseq_mark_stage_inventory.tsv")
    chip_evidence = read_tsv(atlas / "data/chipseq_gene_mark_stage_evidence.tsv")
    chip_positions = read_tsv(atlas / "data/chipseq_gene_position_mark.tsv")
    chip_loci = read_tsv(atlas / "data/chipseq_gene_peak_loci.tsv")
    chip_segments = read_tsv(atlas / "data/chipseq_all_consensus_segments.tsv")
    compatibility = read_tsv(atlas / "data/rna_chip_stage_compatibility.tsv")
    if sum(int(row["n_samples"]) for row in chip_inventory) != manifest["chipseq_samples"]:
        errors.append("ChIP IP library inventory differs from manifest")
    if sum(int(row["n_controls"]) for row in chip_inventory) != manifest["chipseq_controls"]:
        errors.append("ChIP Input library inventory differs from manifest")
    if len({row["stage_or_condition"] for row in chip_inventory}) != manifest["chipseq_stages"]:
        errors.append("ChIP stage count differs from manifest")
    if not chip_evidence or any(row["gene_id"] not in reference_ids for row in chip_evidence):
        errors.append("ChIP evidence is empty or contains genes outside the frozen reference")
    allowed_position_classes = {"promoter_TSS", "gene_body", "distal_intergenic", "other_annotated"}
    if not chip_positions or {row["position_class"] for row in chip_positions} - allowed_position_classes:
        errors.append("all-gene ChIP position evidence is empty or contains an invalid class")
    if any(int(row["n_unique_peaks"]) < 1 for row in chip_positions):
        errors.append("ChIP position evidence contains a non-positive segment count")
    locus_keys = {
        (row["study_id"], row["context_id"], row["gene_id"], row["peak_id"])
        for row in chip_loci
    }
    if not chip_loci or len(locus_keys) != len(chip_loci):
        errors.append("all-gene ChIP locus evidence is empty or contains duplicate segment identities")
    if any(
        row["gene_id"] not in reference_ids
        or int(row["peak_start"]) < 0
        or int(row["peak_end"]) <= int(row["peak_start"])
        or row["position_class"] not in allowed_position_classes
        for row in chip_loci
    ):
        errors.append("ChIP locus evidence contains an invalid gene, coordinate, or position class")
    chip_config = json.loads((root / "analyses/atlas/atlas_config.json").read_text(encoding="utf-8"))["chipseq_integration"]
    expected_native_associations, expected_native_segments = set(), {}
    for entry in chip_config["studies"]:
        source_manifest = root / entry["source_manifest"]
        source = json.loads(source_manifest.read_text(encoding="utf-8"))
        if source["commit"] != entry["commit"]:
            errors.append("ChIP source commit differs from configured frozen version")
        for record in source["files"]:
            if sha256(source_manifest.parent / record["path"]) != record["sha256"]:
                errors.append(f"frozen ChIP source checksum mismatch: {record['path']}")
        result = source_manifest.parent / f"results/chipseq/{entry['id']}"
        regional = result / "manifests/regional_publication.json"
        is_regional = regional.is_file()
        native_directories = [result / "broad/chipseq/peak_annotation" / g["condition"] for g in json.loads(regional.read_text())["groups"]] if is_regional else [result / "broad/chipseq/peak_annotation/peak_annotation_aggregate"]
        if entry.get("admission_scope") == "published_descriptive_catalog_only":
            native_directories += sorted((result / "narrow/chipseq/peak_annotation").glob("*"))
        for native in native_directories:
            suffix = ".gz" if is_regional else ""
            for row in read_tsv(native / ("peak_gene_associations.tsv" + suffix)):
                expected_native_associations.add((entry["id"], row["gene_id"], row["peak_id"], row["category"]))
            for row in read_tsv(native / ("annotated_peaks.tsv" + suffix)):
                expected_native_segments[(entry["id"], row["peak_id"])] = (row["chrom"], row["start"], row["end"], row["category"])
                if is_regional and (row["gene_assignment_ambiguous"] not in {"true","false","not_assessed"} or row["occupancy_status"] != "not_assessed"):
                    errors.append("published regional assessment state was changed")
        study_page = atlas / f"chipseq_{entry['id']}.html"
        if not study_page.is_file() or "atlas_payload.js" in study_page.read_text(encoding="utf-8"):
            errors.append("ChIP study page is missing or unnecessarily loads RNA matrices")
    observed_native_associations = {(r["study_id"], r["gene_id"], r["peak_id"], r["peak_location"]) for r in chip_loci}
    if observed_native_associations != expected_native_associations:
        errors.append("full native gene–segment associations were changed or restricted to selected genes")
    observed_native_segments = {(r["study_id"], r["peak_id"]): (r["chrom"], r["start"], r["end"], r["category"]) for r in chip_segments}
    if observed_native_segments != expected_native_segments:
        errors.append("native consensus segments/intergenic intervals/coordinates were not preserved")
    explorer_path = atlas / "data/chipseq/explorer.json"
    explorer = json.loads(explorer_path.read_text(encoding="utf-8"))
    explorer_segments = [dict(zip(explorer["region_fields"], row)) for row in explorer["consensus"]]
    lazy_cores=[];lazy_support_rows=0;lazy_idr=[]
    for descriptor in explorer.get("lazy_catalogs",[]):
        packed=json.loads((atlas/descriptor["url"]).read_text());c=descriptor["context"]
        for values in packed["rows"]:
            row={**packed["defaults"],**dict(zip(packed["fields"],values))};row["id"]=packed["id_prefix"]+row["id"];explorer_segments.append(row)
        lazy_cores.extend((c["study"],c["stage"],c["factor"],*r) for r in packed["cores"])
        lazy_idr.extend(packed["idr"])
        for path in (atlas/descriptor["support_base"]).parent.glob(Path(descriptor["support_base"]).name+"*.json"):
            lazy_support_rows+=sum(len(rows) for rows in json.loads(path.read_text())["regions"].values())
    if {(r["study"], r["id"]): (r["chrom"], str(r["start"]), str(r["end"]), r["category"]) for r in explorer_segments} != expected_native_segments:
        errors.append("ChIP explorer changed consensus/intergenic source intervals")
    native_peaks = set()
    for entry in chip_config["studies"]:
        frozen = (root / entry["source_manifest"]).parent
        for source in (frozen / f"results/chipseq/{entry['id']}/broad/080-peak-calling").glob("*/peaks.broadPeak"):
            peak_manifest = json.loads((source.parent / "manifest.json").read_text())
            for line in source.read_text().splitlines():
                p = line.split("\t")
                if p and p[0]:native_peaks.add((entry["id"], peak_manifest["record_id"], p[3], p[0], int(p[1]), int(p[2])))
        for source in (frozen / f"results/chipseq/{entry['id']}/broad/chipseq/consensus").rglob("replicate_evidence.tsv"):
            consensus_manifest = json.loads((source.parent / "manifest.json").read_text(encoding="utf-8"))
            reps = {r["evidence_replicate_id"]: r["record_id"] for r in consensus_manifest["replicates"]}
            for row in read_tsv(source):
                native_peaks.add((entry["id"], reps[row["replicate_id"]], row["original_peak_name"], row["chrom"], int(row["start"]), int(row["end"])))
    explorer_peaks = [dict(zip(explorer["region_fields"], row)) for row in explorer["library_peaks"]]
    if {(r["study"], r["library"], r["id"], r["chrom"], r["start"], r["end"]) for r in explorer_peaks} != native_peaks or len(explorer_peaks) != len(native_peaks):
        errors.append("ChIP explorer per-library calls differ from native source")
    if explorer.get("core_segments"):
        expected_cores = set()
        for entry in chip_config["studies"]:
            result = (root / entry["source_manifest"]).parent / f"results/chipseq/{entry['id']}"
            groups=json.loads((result / "manifests/regional_publication.json").read_text())["groups"]
            if entry.get("admission_scope") == "published_descriptive_catalog_only":
                groups+=json.loads((result / "narrow/manifests/native_execution_summary.json").read_text())["groups"]
            for stage in groups:
                with gzip.open(result / f"{'narrow' if stage['target']=='H3K4me3' else 'broad'}/chipseq/consensus/{stage['condition']}/support_core.bed.gz", "rt") as handle:
                    for line in handle:
                        chrom,start,end,*_ = line.rstrip().split("\t")
                        expected_cores.add((entry["id"],stage["condition"],stage["target"],chrom,int(start),int(end)))
        if {tuple(row) for row in explorer["core_segments"]} | set(lazy_cores) != expected_cores:
            errors.append("supplementary shared core differs from native source")
        support = json.loads((atlas / explorer["evidence_url"]).read_text())
        if sum(map(len,support["regions"].values())) != 8049:
            errors.append("regional per-library coverage lost zero-support rows")
        if lazy_support_rows and lazy_support_rows != 484210:errors.append("lazy Picard support lost zero rows")
        if lazy_idr and len(lazy_idr)!=18179:errors.append("native Picard IDR intervals changed")
    if explorer["signal_tracks"] or explorer["signal_status"] != "NOT_INCLUDED_IN_PUBLIC_SNAPSHOT":
        errors.append("ChIP explorer misrepresents unavailable continuous signal")
    for name in ("chipseq.html", "chipseq.js", "chipseq_logic.js", "chipseq_tracks.js", "vendor/chipseq/bbi_reader.js"):
        if not (atlas / name).is_file():
            errors.append(f"Missing ChIP explorer asset: {name}")
    if "atlas_payload.js" in (atlas / "chipseq.html").read_text(encoding="utf-8"):
        errors.append("ChIP explorer unnecessarily loads RNA matrices")
    if any(int(r["support"]) < next(int(c["min_libraries"]) for c in chip_inventory if c["context_id"] == r["context_id"]) for r in chip_segments):
        errors.append("unsupported ChIP consensus segment")
    for row in chip_loci:
        native = expected_native_segments[(row["study_id"], row["peak_id"])]
        if (row["peak_chrom"], row["peak_start"], row["peak_end"], row["peak_location"]) != native:
            errors.append("gene-locus coordinates/classes differ from the native annotation")
            break
    if sum(int(row["n_peaks"]) for row in chip_evidence) != len(chip_loci):
        errors.append("ChIP summary segment counts differ from gene-locus records")
    for legacy_name in ("chipseq_candidate_peak_loci.tsv", "chipseq_candidate_position_mark.tsv"):
        if (atlas / "data" / legacy_name).exists():
            errors.append("retired candidate-only ChIP export remains in the atlas")
    allowed_match_states = {"DIRECT_STAGE_CONTEXT", "CONTEXTUAL_STAGE_FAMILY", "NOT_COVERED"}
    if not compatibility or {row["match_status"] for row in compatibility} - allowed_match_states:
        errors.append("RNA–ChIP stage compatibility uses an invalid state")
    covered = {row["study"] for row in compatibility if row["match_status"] != "NOT_COVERED"}
    if manifest["rna_studies_with_stage_context"] != len(covered):
        errors.append("RNA study stage-context count differs from actual ChIP coverage")

    payload_path = atlas / "data/atlas_payload.js"
    payload_text = payload_path.read_text(encoding="utf-8")
    if not payload_text.startswith("window.HF_ATLAS=") or not payload_text.endswith(";\n"):
        errors.append("browser payload wrapper is invalid")
    else:
        payload = json.loads(payload_text[len("window.HF_ATLAS="):-2])
        if payload.get("locus_display") != display:
            errors.append("browser locus reference windows differ from the recorded display policy")
        if len(payload.get("genes", [])) != EXPECTED["genes"]:
            errors.append("browser payload gene count mismatch")
        if len(payload.get("samples", [])) != EXPECTED["biological_samples"]:
            errors.append("browser payload sample count mismatch")
        if len(payload.get("contrasts", [])) != EXPECTED["contrasts"]:
            errors.append("browser payload contrast count mismatch")
        profiles = {p["study"]: p for p in payload.get("study_profiles", [])}
        for study in ("PRJNA312093", "PRJNA1177852", "PRJNA428470"):
            if not profiles.get(study, {}).get("interpretation_warning"):
                errors.append(f"browser admission warning missing: {study}")
        source_tpm = {row["gene_id"]: row for row in read_tsv(root / "results/PRJNA312093/expression/tpm_matrix.tsv")}
        sample_fields = payload["sample_fields"]
        refit_columns = []
        for column, record in enumerate(payload["samples"]):
            sample = dict(zip(sample_fields, record))
            if sample["study"] == "PRJNA312093":
                refit_columns.append((column, "PRJNA312093__" + sample["sample_id"]))
        for gene_index, gene in enumerate(payload["genes"]):
            if any(payload["tpm"][gene_index][column] != float(f"{float(source_tpm[gene][key]):.5g}") for column, key in refit_columns):
                errors.append(f"browser refit expression differs from source: {gene}")
                break
        if len(payload.get("gene_annotations", [])) != EXPECTED["genes"]:
            errors.append("browser payload annotation count mismatch")
        fields = ["transcript_id", "chromosome", "start", "end", "strand", "exons"]
        expected_browser_models = [[[model[field] for field in fields] for model in expected_models[gene]] for gene in payload["genes"]]
        if payload.get("gene_model_fields") != fields or payload.get("gene_models") != expected_browser_models:
            errors.append("browser exon structures are incomplete or attached to the wrong genes")
        if not {"annotation_name", "description", "biotype", "description_source", "source_accession", "previous_stable_id", "chromosome", "start", "end", "strand"}.issubset(payload.get("gene_annotation_fields", [])):
            errors.append("browser payload lacks WBPS19 gene annotation fields")
        if len(payload.get("chip_mark_stages", [])) != len(chip_inventory):
            errors.append("browser payload ChIP inventory count mismatch")
        delayed = {c["context"]["study"] for c in explorer.get("lazy_catalogs", [])}
        if len(payload.get("chip_evidence", [])) != sum(row["study_id"] not in delayed for row in chip_evidence):
            errors.append("browser eager ChIP evidence count mismatch")
        if len(payload.get("chip_loci", [])) != sum(row["study_id"] not in delayed for row in chip_loci):
            errors.append("browser eager ChIP locus count mismatch")
        # ChIP retains reference genes without RNA measurements. Lazy delivery must
        # preserve all assignments rather than restricting them to the RNA universe.
        lazy_loci = set()
        store = payload.get("chip_gene_store", {})
        for bucket in range(store.get("buckets", 0)):
            packed = json.loads((atlas / store["pattern"].replace("{bucket}", str(bucket))).read_text())
            for values in packed["rows"]:
                row = {**packed["defaults"], **dict(zip(packed["fields"], values))}
                prefix, suffix = row["peak_id"]
                peak = packed["peak_id_prefixes"][prefix] + "Interval_" + suffix
                lazy_loci.add((row["study_id"], row["gene_id"], peak, row["peak_location"]))
        if lazy_loci != {row for row in expected_native_associations if row[0] in delayed}:
            errors.append("lazy ChIP gene assignments differ from the complete native reference catalog")
        expected_summaries = {}
        metadata = [dict(zip(payload["chip_mark_stage_fields"], row)) for row in payload["chip_mark_stages"]]
        context_index = {row["context_id"]: i for i, row in enumerate(metadata)}
        for row in chip_loci:
            if row["study_id"] not in delayed:
                continue
            key = (row["gene_id"], context_index[row["context_id"]], int(row["support"]))
            values = expected_summaries.setdefault(key, [0, 0, 0, 0, 0])
            values[0] += 1
            values[1] += row["promoter_flag"] == "true"
            values[2] += row["gene_assignment_ambiguous"] == "true"
            values[3] += row["gene_assignment_ambiguous"] == "not_assessed"
            values[4] |= 1 << payload["chip_category_values"].index(row["peak_location"])
        observed_summaries = {tuple(row[:3]): row[3:] for row in payload.get("chip_summaries", [])}
        if observed_summaries != expected_summaries:
            errors.append("browser compact ChIP summaries differ from full native assignments")
        if payload.get("chip_scope") != "all_atlas_genes":
            errors.append("browser ChIP scope is not all atlas genes")
        size_report = json.loads((atlas / "size_report.json").read_text(encoding="utf-8"))
        if size_report["payload_bytes"] != payload_path.stat().st_size:
            errors.append("atlas size balance does not match actual browser payload")
        if len(payload.get("external_annotation", [])) != manifest.get("external_annotation_mapped"):
            errors.append("browser payload Anderson evidence count mismatch")
        if len(payload.get("cell_markers", [])) != manifest.get("cell_marker_mapped"):
            errors.append("browser payload Wang evidence count mismatch")
        if payload.get("cell_marker_source_count") != len(markers):
            errors.append("browser Wang source coverage count mismatch")
        if payload.get("cell_marker_mapped_gene_count") != len({row["current_gene_id"] for row in markers if row["mapping_status"] == "MAPPED"}):
            errors.append("browser Wang mapped-gene count mismatch")
        if payload.get("lncrna_source_count") != len(lncrnas):
            errors.append("browser Vasconcelos source count mismatch")
        if len(payload.get("lncrna_gene_links", [])) != manifest.get("vasconcelos_current_gene_links"):
            errors.append("browser Vasconcelos gene-link count mismatch")
        if payload.get("almeida_source_count") != len(legacy):
            errors.append("browser Almeida source count mismatch")
        if len(payload.get("almeida_evidence", [])) != manifest.get("almeida_mapped_records"):
            errors.append("browser Almeida mapped-evidence count mismatch")
        for key in ("lncrna_gene_links", "almeida_evidence"):
            if any(not isinstance(row[0], int) or row[0] < 0 or row[0] >= EXPECTED["genes"] for row in payload.get(key, [])):
                errors.append(f"browser {key} contains an invalid gene index")
        marker_fields = payload.get("cell_marker_fields", [])
        if {"gene_index", "source_record_id", "reported_statistic"}.issubset(marker_fields):
            field_pos = {field: marker_fields.index(field) for field in ("gene_index", "source_record_id", "reported_statistic")}
            gene_pos = {gene: index for index, gene in enumerate(payload.get("genes", []))}
            expected_wang = sorted((gene_pos[row["current_gene_id"]], row["source_record_id"], row["reported_statistic"])
                                   for row in markers if row["mapping_status"] == "MAPPED" and row["current_gene_id"] in gene_pos)
            actual_wang = sorted((row[field_pos["gene_index"]], row[field_pos["source_record_id"]], row[field_pos["reported_statistic"]])
                                 for row in payload.get("cell_markers", []))
            if actual_wang != expected_wang:
                errors.append("browser Wang records are attached to the wrong genes or statistics")
        else:
            errors.append("browser Wang fields are incomplete")

    for svg in sorted((atlas / "figures").glob("*.svg")):
        try:
            root_element = ET.parse(svg).getroot()
            if not root_element.tag.endswith("svg"):
                errors.append(f"not an SVG document: {svg.name}")
        except ET.ParseError as exc:
            errors.append(f"invalid SVG {svg.name}: {exc}")

    report = (atlas / "atlas.html").read_text(encoding="utf-8")
    if "th{position:sticky;top:49px" in report:
        errors.append("table headers can cover a single Wang evidence row")
    required_sections = ('id="overview"', 'id="samples"', 'id="biology"', 'id="integration"', 'id="literature"', 'id="genes"', 'id="data"')
    for marker in required_sections:
        if marker not in report:
            errors.append(f"missing report section: {marker}")
    for marker in ("Gene locus and linked peaks", "External publication evidence", "Published external evidence for group members", "SM_V10 / WBPS19 gene annotation", "Vasconcelos et al.", "Almeida et al."):
        if marker not in report:
            errors.append(f"missing gene-level integration component: {marker}")
    # JavaScript templates contain dynamic links; validate rendered/static HTML here.
    static_report = re.sub(r"<script\b[^>]*>[\s\S]*?</script>", "", report, flags=re.I)
    for link in re.findall(r'(?:src|href)="([^"#]+)"', static_report):
        if "://" in link or link.startswith("mailto:"):
            continue
        local_path = atlas / link.split("?", 1)[0].split("#", 1)[0]
        if not local_path.is_file():
            errors.append(f"broken local report link: {link}")
        elif "?v=" in link and link.rsplit("?v=", 1)[1] != sha256(local_path)[:12]:
            errors.append(f"stale scientific asset version: {link}")

    for path in atlas.rglob("*"):
        if not path.is_file() or path.suffix.lower() in {".png", ".gz"}:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        for token in FORBIDDEN:
            if token in text:
                errors.append(f"private runtime token in {path.relative_to(atlas)}: {token}")
    return errors


def main() -> int:
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[2]
    atlas = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else None
    errors = validate(root, atlas)
    print(json.dumps({"status": "PASS" if not errors else "FAIL", "errors": errors}, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
