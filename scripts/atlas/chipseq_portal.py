"""Export a standalone ChIP explorer without importing the RNA matrices."""
import csv
import hashlib
import json
from pathlib import Path
import shutil
from chipseq_lazy import export_lazy


def build_chipseq_portal(root, output, chip_data, annotations, gene_models, rna_genes):
    reference_path = root / "config/reference/chipseq_chromosomes_wbps19.json"
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    if any(p["reference_id"] != reference["reference_id"] for p in chip_data["profiles"]):
        raise ValueError("ChIP explorer reference mismatch")
    chroms = reference["chromosomes"]
    gene_fields = ["id", "description", "chrom", "start", "end", "strand", "previous", "rna"]
    genes, models = [], {chrom: {} for chrom in chroms}
    for gene, annotation in sorted(annotations.items()):
        chrom, start, end = annotation["chromosome"], int(annotation["start"]) - 1, int(annotation["end"])
        if chrom not in chroms or not 0 <= start < end <= chroms[chrom]:
            raise ValueError("Invalid gene coordinates for ChIP explorer")
        genes.append([gene, annotation.get("functional_name") or annotation.get("description", ""),
                      chrom, start, end, annotation["strand"], annotation.get("previous_stable_id", ""), gene in rna_genes])
        models[chrom][gene] = [[m["transcript_id"], [[s - 1, e] for s, e in m["exons"]]] for m in gene_models.get(gene, [])]
    fields = ["study", "stage", "factor", "library", "id", "chrom", "start", "end", "category", "genes", "support", "signal", "p", "q"]
    segments = [[r["study_id"], r["stage_or_condition"], r["mark_or_factor"], "", r["peak_id"],
                 r["chrom"], int(r["start"]), int(r["end"]), r["category"], r["gene_ids"], r["support"], None, None, None]
                for r in chip_data["segments"]]
    peaks = [[r["study_id"], r["stage_or_condition"], r["mark_or_factor"], r["run_accession"], r["peak_id"],
              r["chrom"], r["start"], r["end"], "not_annotated", "", None, r["signal_value"], r["p_value"], r["q_value"]]
             for r in chip_data["library_peaks"]]
    extra_fields = ["method", "annotation_provider", "assignment_method", "distance_to_tss", "libraries", "total_libraries", "support_category", "core_bp", "ambiguity", "ambiguity_reason", "ambiguity_scope", "alternatives", "occupancy"]
    for row, native in zip(segments, chip_data["segments"]):
        row.extend([native.get("consensus_method","replicate_support"),native.get("annotation_provider","python_interval_v1"),
                    native.get("association_method","overlap_priority_first_gene"),native.get("distance_to_tss",""),
                    native.get("calling_libraries",""),native.get("total_library_count",None),native.get("library_support_category",""),
                    native.get("support_core_bp",None),native.get("gene_assignment_ambiguous","not_assessed"),
                    native.get("gene_assignment_ambiguity_reason",""),native.get("gene_assignment_ambiguity_scope",""),
                    native.get("alternative_gene_ids","[]"),native.get("occupancy_status","not_assessed")])
    for row in peaks:row.extend([None]*len(extra_fields))
    fields += extra_fields
    from hf_chipseq_evidence import LOCUS_FIELDS
    lazy = export_lazy(output, chip_data, fields, segments, LOCUS_FIELDS)
    chip_data['lazy'] = lazy
    delayed = set(lazy['delayed_studies'])
    identities = set()
    for kind, records in (("consensus", segments), ("library", peaks)):
        for row in records:
            if row[5] not in chroms or not 0 <= row[6] < row[7] <= chroms[row[5]]:
                raise ValueError("Invalid native ChIP browser interval")
            key = kind, row[0], row[3], row[4]
            if key in identities:
                raise ValueError("Duplicate ChIP browser interval")
            identities.add(key)
    payload = {"schema_version": "1.0", "reference": reference, "coordinate_system": "0-based-half-open",
               "gene_fields": gene_fields, "genes": genes, "region_fields": fields,
               "consensus": [r for r in segments if r[0] not in delayed], "library_peaks": peaks,
               "lazy_catalogs": lazy['lazy_catalogs'], "summary_fields": lazy['summary_fields'], "gene_summaries": lazy['gene_summaries'], "category_values": lazy['category_values'],
               "total_regions": len(segments),
               "contexts": lazy['contexts'],
               "core_fields": ["study","stage","factor","chrom","start","end"],
               "core_segments": [[r[k] for k in ("study_id","stage_or_condition","mark_or_factor","chrom","start","end")] for r in chip_data.get("core_segments",[]) if r['study_id'] not in delayed],
               "evidence_url": "data/chipseq/library_support.json",
               "profiles": [{k: v for k, v in p.items() if k not in {"source_files","publication"}} for p in chip_data["profiles"]],
               "signal_tracks": [], "signal_status": "NOT_INCLUDED_IN_PUBLIC_SNAPSHOT",
               "pooled_regions_status": "COORDINATES_NOT_INCLUDED_IN_PUBLIC_SNAPSHOT"}
    data = output / "data/chipseq"
    data.mkdir(parents=True, exist_ok=True)
    (data / "explorer.json").write_text(json.dumps(payload, separators=(",", ":"), ensure_ascii=False) + "\n", encoding="utf-8")
    support_fields = ["record_id","sample_id","called_peak_count","covered_bp","covered_fraction","region_bp","region_core_bp"]
    support_rows = {}
    for r in chip_data.get("regional_support",[]):
        if any(r['peak_id'].startswith(study+'.') for study in delayed):continue
        support_rows.setdefault(r["peak_id"],[]).append([r[k] for k in support_fields])
    (data / "library_support.json").write_text(json.dumps({"fields":support_fields,"regions":support_rows},separators=(",", ":")) + "\n",encoding="utf-8")
    for chrom, content in models.items():
        (data / f"models_{chrom}.json").write_text(json.dumps(content, separators=(",", ":")) + "\n", encoding="utf-8")
    with (data / "library_peaks.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(peaks)
    sources = [reference_path, root / "scripts/atlas/chipseq_portal.py"]
    for name in ("chipseq.html", "chipseq.css", "chipseq.js", "chipseq_logic.js", "chipseq_tracks.js", "chipseq_signal_source.js"):
        source = root / "analyses/atlas/chipseq_portal" / name
        sources.append(source)
        if name != "chipseq_signal_source.js":
            shutil.copy2(source, output / name)
    vendor = root / "analyses/atlas/chipseq_portal/vendor"
    shutil.copytree(vendor, output / "vendor/chipseq", dirs_exist_ok=True)
    sources.extend(p for p in vendor.rglob("*") if p.is_file())
    page = output / "chipseq.html"
    text = page.read_text(encoding="utf-8")
    for name in ("chipseq.css", "chipseq_logic.js", "chipseq_tracks.js", "chipseq.js"):
        digest = hashlib.sha256((output / name).read_bytes()).hexdigest()[:12]
        text = text.replace(f'"{name}"', f'"{name}?v={digest}"')
    page.write_text(text, encoding="utf-8")
    return sources
