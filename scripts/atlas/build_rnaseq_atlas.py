#!/usr/bin/env python3
"""Build the descriptive HelixForge-SMansoni RNA-seq atlas.

The atlas combines accepted per-study artifacts for navigation and descriptive
visualization. It never refits differential-expression models and does not
perform batch correction.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
from hf_chipseq_evidence import load_hf_chipseq, METADATA_FIELDS, EVIDENCE_FIELDS, LOCUS_FIELDS, POSITION_FIELDS, render_study_page
from chipseq_portal import build_chipseq_portal
from gene_collections import build_gene_collections
from atlas_design import apply_atlas_design
from external_rna_evidence import ANNOTATION_FIELDS, MARKER_FIELDS, curate
from literature_rna_evidence import (
    LEGACY_FIELDS,
    LNCRNA_FIELDS,
    curate as curate_literature,
)


STUDY_COLORS = {
    "PRJNA602528": "#0f766e",
    "PRJNA597909": "#2563eb",
    "PRJEB14695": "#7c3aed",
    "PRJEB32839": "#db2777",
    "E-MTAB-451": "#ea580c",
    "PRJEB3190": "#65a30d",
    "E-ERAD-478": "#0891b2",
    "PRJNA704646": "#a16207",
    "PRJNA1150591": "#475569",
    "PRJNA312093": "#be123c",
    "PRJNA1177852": "#0369a1",
    "PRJNA428470": "#854d0e",
}
QC_COLORS = {"PASS": "#16a34a", "REVIEW": "#d97706", "FAIL": "#dc2626", "NOT_RECORDED": "#64748b"}
LIFE_CYCLE_STAGES = ("egg", "miracidium", "sporocyst", "cercaria", "schistosomulum", "juvenile", "adult")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_float(value: object) -> float:
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return math.nan


def safe_int(value: object) -> int:
    number = safe_float(value)
    return int(number) if math.isfinite(number) else 0


def compact_number(value: float, digits: int = 5) -> float | None:
    if not math.isfinite(value):
        return None
    return float(f"{value:.{digits}g}")


def stage_family(value: object) -> str:
    """Map life-cycle labels to a conservative cross-assay stage family."""
    text = str(value or "").strip().lower().replace("-", "_").replace(" ", "_")
    aliases = (
        ("miracid", "miracidium"),
        ("sporocyst", "sporocyst"),
        ("cercar", "cercaria"),
        ("schistosom", "schistosomulum"),
        ("juvenile", "juvenile"),
        ("adult", "adult"),
        ("egg", "egg"),
    )
    for token, family in aliases:
        if token in text:
            return family
    return text


def life_cycle_sort_key(sample: dict[str, object]) -> tuple[int, float, float]:
    """Order known stages biologically, then developmental age and recorded time."""
    stage = str(sample.get("stage") or "").strip().lower().replace("-", "_").replace(" ", "_")
    family = stage_family(stage)
    rank = LIFE_CYCLE_STAGES.index(family) if family in LIFE_CYCLE_STAGES else len(LIFE_CYCLE_STAGES)
    age_match = re.search(r"(?:^|_)(?:d(\d+(?:\.\d+)?)|(\d+(?:\.\d+)?)(h|d|wk))(?:_|$)", stage)
    age_hours = math.inf
    if age_match:
        if age_match.group(1):
            age_hours = float(age_match.group(1)) * 24
        else:
            age_hours = float(age_match.group(2)) * {"h": 1, "d": 24, "wk": 168}[age_match.group(3)]
    time_hours = safe_float(sample.get("time_hours"))
    return rank, age_hours, time_hours if math.isfinite(time_hours) else math.inf


def stage_match_status(rna_stage: object, chip_stage: object, rna_family_override: str | None = None) -> str:
    """Classify stage compatibility without implying paired experiments."""
    rna_text = str(rna_stage or "").strip().lower().replace("-", "_").replace(" ", "_")
    rna_family = rna_family_override or stage_family(rna_text)
    chip_family = stage_family(chip_stage)
    if not rna_family or rna_family != chip_family:
        return "NOT_COVERED"
    exact_labels = {
        "adult": {"adult"},
        "cercaria": {"cercaria", "cercariae"},
        "miracidium": {"miracidium", "miracidia"},
        "schistosomulum": {"schistosomulum", "schistosomula"},
        "sporocyst": {"sporocyst", "sporocysts"},
    }
    if rna_family_override is None and rna_text in exact_labels.get(rna_family, set()):
        return "DIRECT_STAGE_CONTEXT"
    return "CONTEXTUAL_STAGE_FAMILY"


def load_matrix(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        genes: list[str] = []
        values: list[list[float]] = []
        for row in reader:
            if not row:
                continue
            genes.append(row[0])
            values.append([safe_float(value) for value in row[1:]])
    matrix = np.asarray(values, dtype=float)
    if matrix.shape != (len(genes), len(header) - 1):
        raise ValueError(f"matrix shape mismatch: {path}")
    if not np.isfinite(matrix).all():
        raise ValueError(f"non-finite expression value: {path}")
    return genes, header[1:], matrix


def select_expression_samples(columns, matrix, allowlist_path, study):
    """Select the frozen sample keys; never infer admission from expression."""
    rows = read_tsv(allowlist_path)
    selected = [row.get("sample_key") or row.get("__rowname", "") for row in rows]
    if not selected or any(not key.startswith(study + "__") for key in selected):
        raise ValueError(f"invalid sample allowlist for {study}")
    if len(set(selected)) != len(selected) or len(set(columns)) != len(columns):
        raise ValueError(f"duplicate sample keys for {study}")
    missing = set(selected) - set(columns)
    if missing:
        raise ValueError(f"approved samples missing for {study}: {sorted(missing)}")
    positions = {key: index for index, key in enumerate(columns)}
    return selected, matrix[:, [positions[key] for key in selected]]


def study_spec(root, entry):
    name = "de_spec.json" if entry["differential_expression"] else "de_spec.blocked.json"
    path = root / "config" / entry["id"] / name
    spec = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    if "contrast_sources" in entry:
        approved = set(entry["contrast_sources"])
        available = {contrast["id"] for contrast in spec.get("contrasts", [])}
        if not approved or approved - available:
            raise ValueError(f"invalid approved contrast selection for {entry['id']}")
        spec["contrasts"] = [c for c in spec["contrasts"] if c["id"] in approved]
        spec["source_analysis_id"] = spec.get("analysis_id", "")
        spec["analysis_id"] = entry["analysis_id"]
    return path, spec


def approved_contrast_sources(root, entry, spec):
    records = entry["contrast_sources"]
    if set(records) != {contrast["id"] for contrast in spec["contrasts"]}:
        raise ValueError(f"approved contrast sources differ from specification for {entry['id']}")
    sources = []
    for contrast in spec["contrasts"]:
        record = records[contrast["id"]]
        path = root / record["path"]
        if not path.is_file() or sha256(path) != record["sha256"]:
            raise ValueError(f"approved contrast checksum mismatch: {record['path']}")
        sources.append((contrast["id"], path))
    return sources


def parse_candidate_genes(path: Path) -> tuple[list[str], dict[str, list[str]]]:
    ordered: list[str] = []
    groups: dict[str, list[str]] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        group, raw_genes = line.split(":", 1)
        genes = [item.strip() for item in raw_genes.split(",") if item.strip()]
        groups[group.strip()] = list(dict.fromkeys(genes))
        for gene in genes:
            if gene not in ordered:
                ordered.append(gene)
    return ordered, groups


def parse_curated_candidate_genes(path: Path) -> tuple[list[str], dict[str, list[str]], dict[str, dict[str, str]], list[dict[str, str]]]:
    """Keep the complete audit table while excluding rejected hypotheses from analysis."""
    records = read_tsv(path)
    required = {"record_id", "input_gene_id", "sm_v10_gene_id", "broad_class", "functional_group", "scope_tier", "curation_status", "mapping_status", "evidence_source", "source_url", "notes"}
    if not records or not required.issubset(records[0]):
        raise ValueError(f"invalid curated candidate table: {path}")
    active_tiers = {"CORE_CURATED", "CORE_ADD", "REVIEW_ADD", "CONTEXT"}
    class_labels = {"WRITER": "Writers", "READER": "Readers", "ERASER": "Erasers", "HP1": "HP1", "CHAPERONE": "Chaperones", "COFACTOR": "Cofactors", "REMODELER": "Remodelers"}
    excluded_gene_ids = {row["input_gene_id"] for row in records if row["scope_tier"] == "EXCLUDED_REVIEW"}
    ordered: list[str] = []
    groups: dict[str, list[str]] = {}
    metadata: dict[str, dict[str, str]] = {}
    seen_records: set[str] = set()
    for row in records:
        if not row["record_id"] or row["record_id"] in seen_records:
            raise ValueError(f"duplicate or missing curated record ID: {row['record_id']}")
        seen_records.add(row["record_id"])
        if row["input_gene_id"] in excluded_gene_ids:
            continue
        if row["scope_tier"] not in active_tiers:
            raise ValueError(f"unknown candidate tier: {row['scope_tier']}")
        gene = row["input_gene_id"]
        resolved = row["sm_v10_gene_id"]
        if not gene or not resolved:
            raise ValueError(f"unmapped active candidate: {row['record_id']}")
        if row["broad_class"] not in class_labels:
            raise ValueError(f"unknown candidate class: {row['broad_class']}")
        group = class_labels[row["broad_class"]]
        groups.setdefault(group, [])
        if gene not in groups[group]:
            groups[group].append(gene)
        if gene not in metadata:
            ordered.append(gene)
            metadata[gene] = {"resolved_gene_id": resolved, "mapping_status": row["mapping_status"]}
        elif metadata[gene]["resolved_gene_id"] != resolved or metadata[gene]["mapping_status"] != row["mapping_status"]:
            raise ValueError(f"conflicting curated mapping for {gene}")
        for key in ("scope_tier", "broad_class", "functional_group", "curation_status", "record_id", "evidence_source", "source_url", "notes"):
            target = {"scope_tier": "scope_tiers", "broad_class": "broad_classes", "functional_group": "functional_groups", "curation_status": "curation_statuses", "record_id": "record_ids", "evidence_source": "evidence_sources", "source_url": "source_urls", "notes": "curation_notes"}[key]
            values = metadata[gene].setdefault(target, "").split(";") if metadata[gene].get(target) else []
            if row[key] and row[key] not in values:
                values.append(row[key])
            metadata[gene][target] = ";".join(values)
    groups = {label: groups[label] for label in class_labels.values() if label in groups}
    return ordered, groups, metadata, records


def load_gene_annotations(path: Path) -> dict[str, dict[str, str]]:
    annotations = {row["gene_id"]: row for row in read_tsv(path)}
    if not annotations:
        raise ValueError(f"gene annotation table is empty: {path}")
    return annotations


def qc_records(path: Path | None, study: str) -> dict[str, dict[str, object]]:
    if path is None:
        return {}
    records: dict[str, dict[str, object]] = {}
    for row in read_tsv(path):
        sample = row["sample_id"]
        classification = row.get("classification", "")
        if not classification:
            classification = "PASS" if row.get("qc_flag") in {"NONE", "PASS"} else row.get("qc_flag", "NOT_RECORDED")
        mapping = row.get("salmon_mapping_percent", "")
        retention = row.get("trim_retention_percent", "")
        records[sample] = {
            "qc_classification": classification or "NOT_RECORDED",
            "trim_retention_percent": safe_float(retention),
            "salmon_mapping_percent": safe_float(mapping),
            "qc_source": path.as_posix(),
        }
    return records


def pca_coordinates(expression: np.ndarray, top_genes: int) -> tuple[np.ndarray, np.ndarray]:
    # expression is genes x samples; PCA observations are biological samples.
    x = np.log2(expression.T + 1.0)
    variances = np.var(x, axis=0)
    eligible = np.flatnonzero(variances > 0)
    if eligible.size < 2:
        raise ValueError("PCA requires at least two variable genes")
    selected = eligible[np.argsort(variances[eligible])[-min(top_genes, eligible.size):]]
    x = x[:, selected]
    x -= np.mean(x, axis=0)
    standard_deviation = np.std(x, axis=0)
    standard_deviation[standard_deviation == 0] = 1.0
    x /= standard_deviation
    u, singular_values, _ = np.linalg.svd(x, full_matrices=False)
    coordinates = u[:, :2] * singular_values[:2]
    eigenvalues = singular_values**2
    explained = eigenvalues[:2] / eigenvalues.sum()
    return coordinates, explained


def xml_text(value: object) -> str:
    return html.escape(str(value), quote=True)


def svg_shell(width: int, height: int, body: str, title: str, description: str) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">'
        f'<title id="title">{xml_text(title)}</title><desc id="desc">{xml_text(description)}</desc>'
        '<rect width="100%" height="100%" fill="#ffffff"/>'
        '<style>text{font-family:Inter,Segoe UI,Arial,sans-serif;fill:#172033}.axis{stroke:#94a3b8;stroke-width:1}'
        '.grid{stroke:#e2e8f0;stroke-width:1}.small{font-size:11px}.label{font-size:13px}.title{font-size:22px;font-weight:700}'
        '.subtitle{font-size:12px;fill:#64748b}</style>'
        f'{body}</svg>'
    )


def scatter_svg(
    rows: list[dict[str, object]],
    color_field: str,
    colors: dict[str, str],
    title: str,
    subtitle: str,
    explained: tuple[float, float],
) -> str:
    width, height = 1100, 760
    left, right, top, bottom = 90, 250, 90, 80
    plot_w, plot_h = width - left - right, height - top - bottom
    xs = [float(row["PC1"]) for row in rows]
    ys = [float(row["PC2"]) for row in rows]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    xpad = (xmax - xmin or 1) * 0.08
    ypad = (ymax - ymin or 1) * 0.08
    xmin, xmax = xmin - xpad, xmax + xpad
    ymin, ymax = ymin - ypad, ymax + ypad

    def sx(value: float) -> float:
        return left + (value - xmin) / (xmax - xmin) * plot_w

    def sy(value: float) -> float:
        return top + plot_h - (value - ymin) / (ymax - ymin) * plot_h

    elements = [
        f'<text class="title" x="{left}" y="38">{xml_text(title)}</text>',
        f'<text class="subtitle" x="{left}" y="60">{xml_text(subtitle)}</text>',
    ]
    for tick in range(6):
        fraction = tick / 5
        x = left + fraction * plot_w
        y = top + fraction * plot_h
        xv = xmin + fraction * (xmax - xmin)
        yv = ymax - fraction * (ymax - ymin)
        elements.extend([
            f'<line class="grid" x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + plot_h}"/>',
            f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}"/>',
            f'<text class="small" text-anchor="middle" x="{x:.1f}" y="{top + plot_h + 22}">{xv:.1f}</text>',
            f'<text class="small" text-anchor="end" x="{left - 10}" y="{y + 4:.1f}">{yv:.1f}</text>',
        ])
    elements.extend([
        f'<line class="axis" x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}"/>',
        f'<text class="label" text-anchor="middle" x="{left + plot_w / 2}" y="{height - 22}">PC1 ({explained[0] * 100:.1f}%)</text>',
        f'<text class="label" text-anchor="middle" transform="translate(24 {top + plot_h / 2}) rotate(-90)">PC2 ({explained[1] * 100:.1f}%)</text>',
    ])
    for row in rows:
        category = str(row[color_field])
        label = f'{row["study"]} | {row["sample_id"]} | {category}'
        elements.append(
            f'<circle cx="{sx(float(row["PC1"])):.2f}" cy="{sy(float(row["PC2"])):.2f}" r="5" '
            f'fill="{colors.get(category, "#64748b")}" fill-opacity="0.82" stroke="#ffffff" stroke-width="1">'
            f'<title>{xml_text(label)}</title></circle>'
        )
    legend_values = list(dict.fromkeys(str(row[color_field]) for row in rows))
    for index, category in enumerate(legend_values):
        y = top + index * 25
        elements.append(f'<circle cx="{left + plot_w + 35}" cy="{y}" r="6" fill="{colors.get(category, "#64748b")}"/>')
        elements.append(f'<text class="small" x="{left + plot_w + 49}" y="{y + 4}">{xml_text(category)}</text>')
    return svg_shell(width, height, "".join(elements), title, subtitle)


def bar_svg(labels: list[str], series: list[tuple[str, list[float], str]], title: str, subtitle: str, y_label: str) -> str:
    width = max(1000, 180 + len(labels) * 70)
    height = 650
    left, right, top, bottom = 90, 40, 100, 170
    plot_w, plot_h = width - left - right, height - top - bottom
    maximum = max((value for _, values, _ in series for value in values), default=1.0) or 1.0
    group_w = plot_w / max(1, len(labels))
    bar_w = min(24, group_w / max(1, len(series)) * 0.75)
    elements = [
        f'<text class="title" x="{left}" y="38">{xml_text(title)}</text>',
        f'<text class="subtitle" x="{left}" y="60">{xml_text(subtitle)}</text>',
    ]
    for tick in range(6):
        value = maximum * tick / 5
        y = top + plot_h - plot_h * tick / 5
        elements.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}"/>')
        elements.append(f'<text class="small" text-anchor="end" x="{left - 8}" y="{y + 4:.1f}">{value:.0f}</text>')
    for index, label in enumerate(labels):
        center = left + (index + 0.5) * group_w
        for series_index, (name, values, color) in enumerate(series):
            value = values[index]
            x = center + (series_index - (len(series) - 1) / 2) * bar_w - bar_w / 2
            bar_h = value / maximum * plot_h
            y = top + plot_h - bar_h
            elements.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w - 1:.1f}" height="{bar_h:.1f}" fill="{color}"><title>{xml_text(name)}: {value:g}</title></rect>')
        elements.append(f'<text class="small" text-anchor="end" transform="translate({center + 4:.1f} {top + plot_h + 14}) rotate(-55)">{xml_text(label)}</text>')
    for index, (name, _, color) in enumerate(series):
        x = left + index * 170
        elements.append(f'<rect x="{x}" y="{height - 35}" width="13" height="13" fill="{color}"/>')
        elements.append(f'<text class="small" x="{x + 19}" y="{height - 24}">{xml_text(name)}</text>')
    elements.append(f'<text class="label" text-anchor="middle" transform="translate(22 {top + plot_h / 2}) rotate(-90)">{xml_text(y_label)}</text>')
    return svg_shell(width, height, "".join(elements), title, subtitle)


def interpolate_color(value: float, low: float, high: float, missing: str = "#e5e7eb") -> str:
    if not math.isfinite(value):
        return missing
    if high <= low:
        return "#ffffff"
    position = max(0.0, min(1.0, (value - low) / (high - low)))
    stops = [(49, 54, 149), (255, 255, 255), (165, 0, 38)]
    if position <= 0.5:
        ratio = position * 2
        a, b = stops[0], stops[1]
    else:
        ratio = (position - 0.5) * 2
        a, b = stops[1], stops[2]
    rgb = tuple(round(a[i] + ratio * (b[i] - a[i])) for i in range(3))
    return "#" + "".join(f"{channel:02x}" for channel in rgb)


def heatmap_svg(
    values: np.ndarray,
    row_labels: list[str],
    column_labels: list[str],
    title: str,
    subtitle: str,
    low: float,
    high: float,
    significance: np.ndarray | None = None,
    marker_label: str = "significant (padj < 0.05 and |log2FC| ≥ 1)",
    row_cluster: tuple[list[int], dict[str, object]] | None = None,
) -> str:
    if row_cluster is not None:
        order, tree = row_cluster
        values = values[order, :]
        row_labels = [row_labels[index] for index in order]
        if significance is not None:
            significance = significance[order, :]
    cell_w = max(14, min(34, 850 // max(1, len(column_labels))))
    cell_h = max(11, min(24, 950 // max(1, len(row_labels))))
    label_width = max((len(label) for label in row_labels), default=20) * 6 + 18
    tree_width = 86 if row_cluster is not None else 0
    left, top, right, bottom = max(220, min(430, label_width)) + tree_width, 125, 80, 250
    width = left + len(column_labels) * cell_w + right
    height = top + len(row_labels) * cell_h + bottom
    elements = [
        f'<text class="title" x="24" y="36">{xml_text(title)}</text>',
        f'<text class="subtitle" x="24" y="58">{xml_text(subtitle)}</text>',
    ]
    if row_cluster is not None and len(row_labels) > 1:
        right_edge = left - 8
        root_height = float(tree["height"])
        positions = {leaf: top + index * cell_h + cell_h / 2 for index, leaf in enumerate(order)}
        branches: list[str] = []

        def draw_branch(node: dict[str, object]) -> tuple[float, float]:
            if "leaf" in node:
                return right_edge, positions[int(node["leaf"])]
            ax, ay = draw_branch(node["left"])
            bx, by = draw_branch(node["right"])
            x = right_edge - (float(node["height"]) / root_height * 70 if root_height > 0 else 0)
            branches.append(f'<path d="M {ax:.1f} {ay:.1f} H {x:.1f} V {by:.1f} H {bx:.1f}"/>')
            return x, (ay + by) / 2

        draw_branch(tree)
        elements.append('<g fill="none" stroke="#47637d" stroke-width="1.2"><title>Average-linkage gene clustering</title>' + "".join(branches) + '</g>')
    for row_index, label in enumerate(row_labels):
        y = top + row_index * cell_h
        elements.append(f'<text class="small" text-anchor="end" x="{left - 8 - tree_width}" y="{y + cell_h * 0.78:.1f}">{xml_text(label)}</text>')
        for column_index in range(len(column_labels)):
            value = float(values[row_index, column_index])
            x = left + column_index * cell_w
            fill = interpolate_color(value, low, high)
            elements.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{fill}" stroke="#ffffff" stroke-width="0.35"><title>{xml_text(label)} | {xml_text(column_labels[column_index])}: {value:.4g}</title></rect>')
            if significance is not None and bool(significance[row_index, column_index]):
                elements.append(f'<circle cx="{x + cell_w / 2:.1f}" cy="{y + cell_h / 2:.1f}" r="{max(1.5, min(cell_w, cell_h) * 0.12):.1f}" fill="#111827"/>')
    for index, label in enumerate(column_labels):
        x = left + index * cell_w + cell_w * 0.62
        y = top + len(row_labels) * cell_h + 10
        elements.append(f'<text class="small" text-anchor="end" transform="translate({x:.1f} {y}) rotate(-58)">{xml_text(label)}</text>')
    legend_x = left
    legend_y = height - 45
    for step in range(101):
        value = low + (high - low) * step / 100
        elements.append(f'<rect x="{legend_x + step * 2}" y="{legend_y}" width="2" height="12" fill="{interpolate_color(value, low, high)}"/>')
    elements.append(f'<text class="small" x="{legend_x}" y="{legend_y + 30}">{low:g}</text>')
    elements.append(f'<text class="small" text-anchor="end" x="{legend_x + 202}" y="{legend_y + 30}">{high:g}</text>')
    if significance is not None:
        elements.append(f'<circle cx="{legend_x + 250}" cy="{legend_y + 6}" r="3" fill="#111827"/><text class="small" x="{legend_x + 260}" y="{legend_y + 10}">{xml_text(marker_label)}</text>')
    return svg_shell(width, height, "".join(elements), title, subtitle)


def average_linkage_cluster(values: np.ndarray) -> tuple[list[int], dict[str, object]]:
    """Deterministic average-linkage clustering of gene expression rows."""
    count = values.shape[0]
    if count == 0:
        return [], {}
    clusters: dict[int, dict[str, object]] = {
        index: {"leaf": index, "size": 1, "min": index, "height": 0.0}
        for index in range(count)
    }
    active = list(range(count))
    key = lambda a, b: (min(a, b), max(a, b))
    distances = {
        (i, j): float(np.sqrt(np.mean((values[i] - values[j]) ** 2)))
        for i in range(count) for j in range(i + 1, count)
    }
    next_id = count
    while len(active) > 1:
        _, _, _, a, b = min(
            (distances[key(a, b)], min(clusters[a]["min"], clusters[b]["min"]),
             max(clusters[a]["min"], clusters[b]["min"]), a, b)
            for i, a in enumerate(active) for b in active[i + 1:]
        )
        first, second = (a, b) if clusters[a]["min"] < clusters[b]["min"] else (b, a)
        merged = {
            "left": clusters[first], "right": clusters[second],
            "size": clusters[a]["size"] + clusters[b]["size"],
            "min": min(clusters[a]["min"], clusters[b]["min"]),
            "height": distances[key(a, b)],
        }
        for other in active:
            if other not in (a, b):
                distances[key(next_id, other)] = (
                    clusters[a]["size"] * distances[key(a, other)]
                    + clusters[b]["size"] * distances[key(b, other)]
                ) / merged["size"]
        clusters[next_id] = merged
        active = [index for index in active if index not in (a, b)] + [next_id]
        next_id += 1
    order: list[int] = []

    def visit(node: dict[str, object]) -> None:
        if "leaf" in node:
            order.append(int(node["leaf"]))
        else:
            visit(node["left"])
            visit(node["right"])

    visit(clusters[active[0]])
    return order, clusters[active[0]]


def categorical_heatmap_svg(
    values: np.ndarray,
    row_labels: list[str],
    column_labels: list[str],
    title: str,
    subtitle: str,
) -> str:
    """Render RNA–ChIP stage coverage without treating categories as quantities."""
    colors = {0: "#e5e7eb", 1: "#f59e0b", 2: "#0f766e"}
    labels = {0: "not covered", 1: "stage-family context", 2: "direct stage context"}
    cell_w = max(42, min(90, 850 // max(1, len(column_labels))))
    cell_h = 25
    left, top, right, bottom = 270, 125, 80, 160
    width = max(760, left + len(column_labels) * cell_w + right)
    height = top + len(row_labels) * cell_h + bottom
    elements = [
        f'<text class="title" x="24" y="36">{xml_text(title)}</text>',
        f'<text class="subtitle" x="24" y="58">{xml_text(subtitle)}</text>',
    ]
    for row_index, label in enumerate(row_labels):
        y = top + row_index * cell_h
        elements.append(f'<text class="small" text-anchor="end" x="{left - 8}" y="{y + 18}">{xml_text(label)}</text>')
        for column_index, column_label in enumerate(column_labels):
            code = int(values[row_index, column_index])
            x = left + column_index * cell_w
            elements.append(
                f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{colors[code]}" '
                f'stroke="#ffffff"><title>{xml_text(label)} | {xml_text(column_label)}: {xml_text(labels[code])}</title></rect>'
            )
    for index, label in enumerate(column_labels):
        x = left + index * cell_w + cell_w * 0.55
        y = top + len(row_labels) * cell_h + 10
        elements.append(f'<text class="small" text-anchor="end" transform="translate({x:.1f} {y}) rotate(-45)">{xml_text(label)}</text>')
    legend_y = height - 40
    cursor = 24
    for code in (2, 1, 0):
        elements.append(f'<rect x="{cursor}" y="{legend_y}" width="14" height="14" rx="2" fill="{colors[code]}"/>')
        elements.append(f'<text class="small" x="{cursor + 20}" y="{legend_y + 12}">{xml_text(labels[code])}</text>')
        cursor += 230
    return svg_shell(width, height, "".join(elements), title, subtitle)


def chip_library_inventory_svg(rows: list[dict[str, object]]) -> str:
    """Show discrete library counts on a shared zero-based axis, without pooling studies."""
    width, left, top, row_h, bottom = 900, 330, 104, 62, 65
    height = top + len(rows) * row_h + bottom
    plot_w = width - left - 80
    maximum = max(1, max(int(row[field]) for row in rows for field in ("n_samples", "n_controls")))
    elements = ['<text class="title" x="24" y="35">HF-processed ChIP-seq libraries</text>',
                '<text class="subtitle" x="24" y="58">Deposited library counts; source independence is reviewed per study.</text>']
    for tick in range(maximum + 1):
        x = left + plot_w * tick / maximum
        elements.append(f'<line class="grid" x1="{x}" x2="{x}" y1="{top-12}" y2="{height-bottom}"/><text class="small" x="{x}" y="{height-bottom+20}" text-anchor="middle">{tick}</text>')
    for index, row in enumerate(rows):
        y = top + index * row_h
        elements.append(f'<text class="label" x="{left-15}" y="{y+12}" text-anchor="end">{xml_text(row["stage_or_condition"])}</text><text class="small" x="{left-15}" y="{y+29}" text-anchor="end">{xml_text(row["study_id"])} · {xml_text(row["mark_or_factor"])}</text>')
        for offset, field, label, color in ((0, "n_samples", "IP", "#0f766e"), (19, "n_controls", "Input", "#7893a7")):
            count = int(row[field])
            bar_w = plot_w * count / maximum
            elements.append(f'<rect x="{left}" y="{y+offset}" width="{bar_w}" height="14" rx="2" fill="{color}"><title>{xml_text(row["study_id"])} · {xml_text(row["stage_or_condition"])} · {label}: {count}</title></rect><text class="small" x="{left+bar_w+8}" y="{y+offset+12}">{label} {count}</text>')
    elements.append(f'<text class="small" x="{left}" y="{height-12}">Number of libraries</text>')
    return svg_shell(width, height, "".join(elements), "HF-processed ChIP-seq libraries", "Native IP and Input counts by study and stage")


def version_report_assets(report: str, output: Path) -> str:
    """Change static asset URLs when their bytes change, preventing stale scientific figures."""
    versions = {}
    def version(match: re.Match) -> str:
        prefix, relative, suffix = match.groups()
        path = output / relative
        if "${" in relative or not path.is_file():
            return match.group(0)
        if relative not in versions:
            versions[relative] = sha256(path)[:12]
        return f'{prefix}{relative}?v={versions[relative]}{suffix}'
    return re.sub(r'((?:src|href)=")([^"?]+\.(?:svg|js|css))(\")', version, report)


def gene_position_mark_map_svg(
    rows: list[dict[str, object]],
    genes: list[str],
    gene_labels: list[str],
) -> str:
    """Facet candidate peak evidence by life stage and genomic position."""
    stages = list(dict.fromkeys(str(row["stage_or_condition"]) for row in rows)) or ["No associated segments"]
    positions = ["promoter_TSS", "gene_body", "distal_intergenic", "other_annotated"]
    marks = list(dict.fromkeys(str(row["mark_or_factor"]) for row in rows)) or ["HP1"]
    mark_colors = {
        "H3K27me3": "#2563eb",
        "H3K4me3": "#16a34a",
        "H3K9ac": "#f59e0b",
        "H3K9me3": "#dc2626",
        "unknown_ChIP": "#0f766e",
        "HP1": "#0f766e",
    }
    left, top, row_h, facet_w, bottom = 340, 105, 25, 245, 205
    width = left + len(stages) * facet_w + 35
    height = top + len(genes) * row_h + bottom
    x_offsets = {position: 34 + index * 58 for index, position in enumerate(positions)}
    mark_jitter = {mark: (index - (len(marks) - 1) / 2) * 4.2 for index, mark in enumerate(marks)}
    elements = [
        '<text class="title" x="24" y="36">Candidate gene–position–mark associations</text>',
        '<text class="subtitle" x="24" y="58">Reaggregated from peak-level evidence; color = mark, size = unique peaks, triangle = promoter/TSS.</text>',
    ]
    for stage_index, stage in enumerate(stages):
        x0 = left + stage_index * facet_w
        elements.append(f'<rect x="{x0}" y="{top}" width="{facet_w - 8}" height="{len(genes) * row_h}" fill="{("#f8fafc" if stage_index % 2 else "#ffffff")}"/>')
        elements.append(f'<text class="label" text-anchor="middle" x="{x0 + (facet_w - 8) / 2}" y="{top - 20}" font-weight="700">{xml_text(stage)}</text>')
        for position in positions:
            x = x0 + x_offsets[position]
            elements.append(f'<line class="grid" x1="{x}" y1="{top}" x2="{x}" y2="{top + len(genes) * row_h}"/>')
            label = position.replace("_", "/")
            label_y = top + len(genes) * row_h + 14
            elements.append(f'<text class="small" text-anchor="end" transform="translate({x + 4} {label_y}) rotate(-48)">{xml_text(label)}</text>')
    for gene_index_value, (gene, label) in enumerate(zip(genes, gene_labels)):
        y = top + gene_index_value * row_h + row_h / 2
        elements.append(f'<line class="grid" x1="{left}" y1="{y + row_h / 2}" x2="{left + len(stages) * facet_w - 8}" y2="{y + row_h / 2}"/>')
        elements.append(f'<text class="small" text-anchor="end" x="{left - 10}" y="{y + 4}">{xml_text(label)}</text>')
    gene_rows = {gene: index for index, gene in enumerate(genes)}
    for row in rows:
        gene = str(row["gene_id"])
        stage = str(row["stage_or_condition"])
        mark = str(row["mark_or_factor"])
        position = str(row["position_class"])
        if gene not in gene_rows or stage not in stages or mark not in marks or position not in positions:
            continue
        count = int(row["n_unique_peaks"])
        x = left + stages.index(stage) * facet_w + x_offsets[position] + mark_jitter[mark]
        y = top + gene_rows[gene] * row_h + row_h / 2
        radius = min(10.5, 2.8 + math.sqrt(max(1, count)) * 1.7)
        title = f"{gene} | {stage} | {mark} | {position}: {count} unique peaks"
        if position == "promoter_TSS":
            points = f"{x:.1f},{y-radius:.1f} {x-radius:.1f},{y+radius:.1f} {x+radius:.1f},{y+radius:.1f}"
            elements.append(f'<polygon points="{points}" fill="{mark_colors[mark]}" fill-opacity="0.82"><title>{xml_text(title)}</title></polygon>')
        else:
            elements.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius:.1f}" fill="{mark_colors[mark]}" fill-opacity="0.82"><title>{xml_text(title)}</title></circle>')
    legend_y = height - 48
    cursor = left
    for mark in marks:
        elements.append(f'<circle cx="{cursor}" cy="{legend_y}" r="5" fill="{mark_colors[mark]}"/>')
        elements.append(f'<text class="small" x="{cursor + 11}" y="{legend_y + 4}">{xml_text(mark)}</text>')
        cursor += 150
    size_x = left + 780
    for index, count in enumerate((1, 5, 20)):
        radius = min(10.5, 2.8 + math.sqrt(count) * 1.7)
        x = size_x + index * 70
        elements.append(f'<circle cx="{x}" cy="{legend_y}" r="{radius:.1f}" fill="#334155" fill-opacity="0.82"/>')
        elements.append(f'<text class="small" x="{x + 14}" y="{legend_y + 4}">{count}</text>')
    elements.append(f'<text class="small" x="{size_x}" y="{legend_y - 22}">Unique peaks</text>')
    return svg_shell(width, height, "".join(elements), "Candidate gene–position–mark associations", "Peak-level candidate ChIP-seq evidence")


def table_html(fields: list[str], rows: list[dict[str, object]], limit: int | None = None) -> str:
    shown = rows if limit is None else rows[:limit]
    head = "".join(f"<th>{html.escape(field.replace('_', ' ').title())}</th>" for field in fields)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(row.get(field, '')))}</td>" for field in fields) + "</tr>"
        for row in shown
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def build_html(template: str, replacements: dict[str, str]) -> str:
    for key, value in replacements.items():
        template = template.replace("{{" + key + "}}", value)
    unresolved = [token for token in template.split("{{")[1:] if "}}" in token]
    if unresolved:
        raise ValueError(f"unresolved HTML placeholders: {unresolved[:3]}")
    return template


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = args.root.resolve()
    config_path = (root / args.config).resolve() if not args.config.is_absolute() else args.config.resolve()
    output = (root / args.output).resolve() if not args.output.is_absolute() else args.output.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if output.exists():
        shutil.rmtree(output)
    data_dir = output / "data"
    figure_dir = output / "figures"
    data_dir.mkdir(parents=True)
    figure_dir.mkdir(parents=True)

    studies = [entry["id"] for entry in config["studies"]]
    annotation_path = root / config["gene_annotations"]
    input_files = [
        config_path,
        root / config["selected_genes"],
        annotation_path,
        Path(__file__).resolve(),
    ]
    all_genes: list[str] | None = None
    sample_keys: list[str] = []
    expression_blocks: list[np.ndarray] = []
    sample_rows: list[dict[str, object]] = []
    study_column_slices: dict[str, slice] = {}
    study_sample_counts: dict[str, int] = {}

    column_start = 0
    for study_entry in config["studies"]:
        study = study_entry["id"]
        tpm_path = root / "results" / study / "expression/tpm_matrix.tsv"
        metadata_path = root / "metadata" / study / "samples.tsv"
        input_files.extend([tpm_path, metadata_path])
        genes, columns, tpm = load_matrix(tpm_path)
        if study_entry.get("sample_allowlist"):
            allowlist_path = root / study_entry["sample_allowlist"]
            input_files.append(allowlist_path)
            columns, tpm = select_expression_samples(columns, tpm, allowlist_path, study)
        if study_entry.get("admission_manifest"):
            input_files.append(root / study_entry["admission_manifest"])
        if all_genes is None:
            all_genes = genes
        elif genes != all_genes:
            raise ValueError(f"gene universe/order differs for {study}")
        metadata = {row["sample_id"]: row for row in read_tsv(metadata_path)}
        qc_path = root / study_entry["qc_table"] if study_entry.get("qc_table") else None
        if qc_path:
            input_files.append(qc_path)
        qc = qc_records(qc_path, study)
        for column in columns:
            prefix = study + "__"
            if not column.startswith(prefix):
                raise ValueError(f"sample column does not use study prefix: {column}")
            sample_id = column[len(prefix):]
            if sample_id not in metadata:
                raise ValueError(f"missing metadata for {column}")
            meta = metadata[sample_id]
            sample_qc = qc.get(sample_id, {})
            row: dict[str, object] = {
                "study": study,
                "sample_id": sample_id,
                "sample_key": column,
                "source_sample_name": meta.get("source_sample_name", "") or meta.get("library_name", "") or meta.get("file_prefix", ""),
                "biosample": meta.get("biosample", ""),
                "condition": meta.get("condition", ""),
                "stage": meta.get("stage", ""),
                "tissue": meta.get("tissue", ""),
                "sex": meta.get("sex", ""),
                "infection_mode": meta.get("infection_mode", ""),
                "treatment": meta.get("treatment", ""),
                "time_hours": meta.get("time_hours", ""),
                "batch": meta.get("batch", "not_declared"),
                "replicate": meta.get("replicate", ""),
                "technical_runs": int(meta.get("technical_runs", "") or meta.get("run_count", "1")),
                "qc_classification": sample_qc.get("qc_classification", "NOT_RECORDED"),
                "trim_retention_percent": compact_number(float(sample_qc.get("trim_retention_percent", math.nan))),
                "salmon_mapping_percent": compact_number(float(sample_qc.get("salmon_mapping_percent", math.nan))),
                "reference_id": config["reference_id"],
            }
            sample_rows.append(row)
        expression_blocks.append(tpm)
        sample_keys.extend(columns)
        study_sample_counts[study] = len(columns)
        study_column_slices[study] = slice(column_start, column_start + len(columns))
        column_start += len(columns)

    assert all_genes is not None
    annotations = load_gene_annotations(annotation_path)
    missing_annotations = [gene for gene in all_genes if gene not in annotations]
    if missing_annotations:
        raise ValueError(
            f"WBPS19 annotation is missing {len(missing_annotations)} atlas genes; "
            f"first={missing_annotations[:5]}"
        )
    expression = np.concatenate(expression_blocks, axis=1)
    gene_index = {gene: index for index, gene in enumerate(all_genes)}
    model_path = root / config["gene_models"]
    reference_path = root / "config/reference/reference_manifest.json"
    input_files.extend([model_path, reference_path, root / "scripts/reference/extract_transcript_models.py"])
    model_source = json.loads(model_path.read_text(encoding="utf-8"))
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    frozen_gff = next(row for row in reference["artifacts"] if row["role"] == "annotation_gff3")
    if (model_source["reference_id"] != config["reference_id"]
            or model_source["coordinate_system"] != "1-based-inclusive"
            or model_source["source_sha256"] != frozen_gff["sha256"]):
        raise ValueError("transcript structures differ from the frozen reference identity")
    gene_models = {gene: model_source["genes"].get(gene, []) for gene in all_genes}
    for gene, models in gene_models.items():
        annotation = annotations[gene]
        if not models and annotation["biotype"] != "pseudogene":
            raise ValueError(f"missing transcript structure for {gene}")
        for model in models:
            if (model["chromosome"] != annotation["chromosome"] or model["strand"] != annotation["strand"]
                    or not int(annotation["start"]) <= model["start"] <= model["end"] <= int(annotation["end"])):
                raise ValueError(f"transcript structure disagrees with gene annotation: {gene}")
    model_export = {**{key: value for key, value in model_source.items() if key != "genes"}, "genes": gene_models}
    (data_dir / "gene_models.json").write_text(json.dumps(model_export, separators=(",", ":")) + "\n", encoding="utf-8")
    previous_id_index: dict[str, list[str]] = defaultdict(list)
    for current_gene, annotation in annotations.items():
        raw_previous = annotation.get("previous_stable_id", "")
        for previous_gene in raw_previous.replace(",", ";").split(";"):
            previous_gene = previous_gene.strip()
            if previous_gene and current_gene in gene_index:
                previous_id_index[previous_gene].append(current_gene)
    sample_fields = [
        "study", "sample_id", "sample_key", "source_sample_name", "biosample", "condition", "stage", "tissue",
        "sex", "infection_mode", "treatment", "time_hours", "batch", "replicate", "technical_runs",
        "qc_classification", "trim_retention_percent", "salmon_mapping_percent", "reference_id",
    ]
    write_tsv(data_dir / "sample_inventory.tsv", sample_fields, sample_rows)

    annotation_fields = [
        "gene_id", "display_label", "functional_name", "annotation_name",
        "description", "biotype", "description_source", "source_accession",
        "previous_stable_id", "chromosome", "start", "end", "strand",
    ]
    annotation_rows = [annotations[gene] for gene in all_genes]
    write_tsv(data_dir / "gene_annotations.tsv", annotation_fields, annotation_rows)

    external_annotation_rows, cell_marker_rows, external_sources = curate(
        root / "analyses/atlas/local_inputs/external_rna", set(all_genes), annotations
    )
    input_files.extend(external_sources)
    input_files.append(root / "scripts/atlas/external_rna_evidence.py")
    write_tsv(data_dir / "external_annotation_evidence.tsv", ANNOTATION_FIELDS, external_annotation_rows)
    write_tsv(data_dir / "cell_type_marker_evidence.tsv", MARKER_FIELDS, cell_marker_rows)

    lncrna_rows, legacy_transcript_rows, literature_sources = curate_literature(
        root / "analyses/atlas/local_inputs/literature_rna", set(all_genes), annotations
    )
    input_files.extend(literature_sources)
    input_files.append(root / "scripts/atlas/literature_rna_evidence.py")
    write_tsv(data_dir / "lncrna_annotation_evidence.tsv", LNCRNA_FIELDS, lncrna_rows)
    write_tsv(data_dir / "legacy_transcript_evidence.tsv", LEGACY_FIELDS, legacy_transcript_rows)

    lncrna_gene_links: list[dict[str, object]] = []
    for row in lncrna_rows:
        for gene in str(row["associated_current_gene_id"]).split(";"):
            if gene:
                lncrna_gene_links.append({**row, "current_gene_id": gene})
    almeida_mapped_rows = [row for row in legacy_transcript_rows if row["current_gene_id"]]

    study_rows: list[dict[str, object]] = []
    for study_entry in config["studies"]:
        study = study_entry["id"]
        selected = [row for row in sample_rows if row["study"] == study]
        qc_counts = Counter(str(row["qc_classification"]) for row in selected)
        mapping_values = [float(row["salmon_mapping_percent"]) for row in selected if row["salmon_mapping_percent"] is not None]
        retention_values = [float(row["trim_retention_percent"]) for row in selected if row["trim_retention_percent"] is not None]
        aggregate_qc_path = root / "results" / study / "qc/qc_summary.json"
        aggregate_qc = {}
        if aggregate_qc_path.is_file() and not study_entry.get("sample_allowlist"):
            input_files.append(aggregate_qc_path)
            aggregate_qc = json.loads(aggregate_qc_path.read_text(encoding="utf-8"))

        def aggregate_metric(values: list[float], suffix: str, operation: str) -> float | None:
            candidates = [
                f"{suffix}_{operation}_percent",
                f"{suffix}_percent_{operation}",
                f"{suffix}_{operation}",
            ]
            for candidate in candidates:
                value = safe_float(aggregate_qc.get(candidate))
                if math.isfinite(value):
                    return compact_number(value)
            if values:
                function = {"min": min, "max": max, "mean": lambda items: sum(items) / len(items)}[operation]
                return compact_number(float(function(values)))
            return None

        study_rows.append({
            "study": study,
            "biological_samples": len(selected),
            "technical_runs": sum(int(row["technical_runs"]) for row in selected),
            "conditions": len({row["condition"] for row in selected}),
            "stages": len({row["stage"] for row in selected}),
            "qc_pass": qc_counts["PASS"],
            "qc_review": qc_counts["REVIEW"],
            "qc_fail": qc_counts["FAIL"],
            "qc_not_recorded": qc_counts["NOT_RECORDED"],
            "mapping_mean_percent": aggregate_metric(mapping_values, "mapping", "mean"),
            "mapping_min_percent": aggregate_metric(mapping_values, "mapping", "min"),
            "mapping_max_percent": aggregate_metric(mapping_values, "mapping", "max"),
            "trim_retention_mean_percent": aggregate_metric(retention_values, "trim_retention", "mean"),
            "trim_retention_min_percent": aggregate_metric(retention_values, "trim_retention", "min"),
            "trim_retention_max_percent": aggregate_metric(retention_values, "trim_retention", "max"),
            "differential_expression": "YES" if study_entry["differential_expression"] else "NO_DESCRIPTIVE_ONLY",
            "interpretation_warning": study_entry.get("interpretation_warning", ""),
        })
    study_fields = [
        "study", "biological_samples", "technical_runs", "conditions", "stages", "qc_pass", "qc_review", "qc_fail",
        "qc_not_recorded", "mapping_mean_percent", "mapping_min_percent", "mapping_max_percent",
        "trim_retention_mean_percent", "trim_retention_min_percent", "trim_retention_max_percent", "differential_expression", "interpretation_warning",
    ]
    write_tsv(data_dir / "study_summary.tsv", study_fields, study_rows)

    # Only reviewed native HF results enter the active atlas; no legacy summaries.
    chip_config = config.get("chipseq_integration")
    chip_metadata_rows, chip_evidence_rows, chip_position_rows, chip_locus_rows = [], [], [], []
    chip_profiles, chip_all_segments, stage_compatibility_rows = [], [], []
    if chip_config:
        chip_data = load_hf_chipseq(root, {**chip_config, "reference_id": config["reference_id"]}, annotations)
        chip_metadata_rows, chip_evidence_rows = chip_data["metadata"], chip_data["evidence"]
        chip_position_rows, chip_locus_rows = chip_data["positions"], chip_data["loci"]
        chip_profiles, chip_all_segments = chip_data["profiles"], chip_data["segments"]
        delayed_studies = {p["study_id"] for p in chip_profiles if p.get("admission_scope") == "published_descriptive_catalog_only"}
        input_files.extend(chip_data["inputs"])
        input_files.extend([root / "scripts/atlas/hf_chipseq_evidence.py", root / "scripts/atlas/hf_chipseq_regional.py", root / "scripts/atlas/hf_chipseq_catalog.py", root / "scripts/atlas/chipseq_lazy.py"])
        input_files.extend(build_chipseq_portal(root, output, chip_data, annotations, gene_models, set(all_genes)))
        for profile in chip_profiles:
            study_id = profile["study_id"]
            (data_dir / f"chipseq_{study_id}.json").write_text(json.dumps(profile, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            write_tsv(data_dir / f"chipseq_{study_id}_qc.tsv", list(profile["samples"][0]), profile["samples"])
            (output / f"chipseq_{study_id}.html").write_text(render_study_page(profile), encoding="utf-8")

        chip_stages = list(dict.fromkeys(str(row["stage_or_condition"]) for row in chip_metadata_rows))
        family_overrides = chip_config.get("rna_stage_family_overrides", {})
        for study in studies:
            rna_stages = list(dict.fromkeys(str(row["stage"]) for row in sample_rows if row["study"] == study))
            for rna_stage in rna_stages:
                for chip_stage in chip_stages:
                    family_override = family_overrides.get(study, {}).get(rna_stage)
                    stage_compatibility_rows.append({
                        "study": study,
                        "rna_stage": rna_stage,
                        "chip_stage": chip_stage,
                        "match_status": stage_match_status(rna_stage, chip_stage, family_override),
                        "interpretation": "independent_experiments_no_sample_pairing",
                    })

        write_tsv(data_dir / "chipseq_mark_stage_inventory.tsv", METADATA_FIELDS, chip_metadata_rows)
        write_tsv(data_dir / "chipseq_gene_mark_stage_evidence.tsv", EVIDENCE_FIELDS, chip_evidence_rows)
        write_tsv(data_dir / "chipseq_gene_position_mark.tsv", POSITION_FIELDS, chip_position_rows)
        write_tsv(data_dir / "chipseq_gene_peak_loci.tsv", LOCUS_FIELDS, chip_locus_rows)
        write_tsv(data_dir / "chipseq_all_consensus_segments.tsv", list(chip_all_segments[0]), chip_all_segments)
        write_tsv(data_dir / "chipseq_regional_candidates.tsv", list(chip_all_segments[0]), chip_all_segments)
        if chip_data.get("idr_segments"):
            write_tsv(data_dir / "chipseq_idr.tsv", list(chip_data["idr_segments"][0]), chip_data["idr_segments"])
        write_tsv(data_dir / "chipseq_regional_library_support.tsv", list(chip_data["regional_support"][0]), chip_data["regional_support"])
        with (data_dir / "chipseq_shared_core.bed").open("w",encoding="utf-8") as core_handle:
            for i,r in enumerate(chip_data["core_segments"]):
                core_handle.write(f'{r["chrom"]}\t{r["start"]}\t{r["end"]}\t{r["study_id"]}|{r["stage_or_condition"]}|{r["mark_or_factor"]}|core_{i+1}\n')

        write_tsv(data_dir / "rna_chip_stage_compatibility.tsv", ["study", "rna_stage", "chip_stage", "match_status", "interpretation"], stage_compatibility_rows)
        # Retire the candidate-only preliminary exports, including stale copies.
        for legacy_name in ("chipseq_candidate_position_mark.tsv", "chipseq_candidate_peak_loci.tsv", "selected_gene_rna_chip_summary.tsv"):
            (data_dir / legacy_name).unlink(missing_ok=True)
        for legacy_name in ("candidate_gene_position_mark_map.svg", "selected_gene_chipseq_evidence.svg"):
            (figure_dir / legacy_name).unlink(missing_ok=True)

        coverage_rows: list[str] = []
        coverage_values: list[list[int]] = []
        for study in studies:
            for rna_stage in dict.fromkeys(str(row["stage"]) for row in sample_rows if row["study"] == study):
                coverage_rows.append(f"{study} · {rna_stage}")
                statuses = {
                    str(row["chip_stage"]): str(row["match_status"])
                    for row in stage_compatibility_rows
                    if row["study"] == study and row["rna_stage"] == rna_stage
                }
                coverage_values.append([
                    2 if statuses[stage] == "DIRECT_STAGE_CONTEXT" else 1 if statuses[stage] == "CONTEXTUAL_STAGE_FAMILY" else 0
                    for stage in chip_stages
                ])
        matched_indices = [index for index, values in enumerate(coverage_values) if any(values)]
        (figure_dir / "rna_chip_stage_coverage.svg").write_text(
            categorical_heatmap_svg(
                np.asarray([coverage_values[index] for index in matched_indices], dtype=int),
                [coverage_rows[index] for index in matched_indices], chip_stages,
                "RNA–ChIP stage coverage",
                "Stage-compatible RNA contexts only; full coverage is exported. Assays are independent.",
            ),
            encoding="utf-8",
        )

        (figure_dir / "chipseq_assay_inventory.svg").write_text(
            chip_library_inventory_svg(chip_metadata_rows),
            encoding="utf-8",
        )

    # PCA: global and within-study views. No adjustment or inference is applied.
    pca_rows: list[dict[str, object]] = []
    global_coordinates, global_explained = pca_coordinates(expression, int(config["pca_top_variable_genes"]))
    global_plot_rows: list[dict[str, object]] = []
    for index, sample in enumerate(sample_rows):
        row = {
            "scope": "GLOBAL_EXPLORATORY",
            "study": sample["study"],
            "sample_id": sample["sample_id"],
            "sample_key": sample["sample_key"],
            "condition": sample["condition"],
            "stage": sample["stage"],
            "PC1": float(global_coordinates[index, 0]),
            "PC2": float(global_coordinates[index, 1]),
            "PC1_variance_percent": float(global_explained[0] * 100),
            "PC2_variance_percent": float(global_explained[1] * 100),
        }
        pca_rows.append(row)
        global_plot_rows.append(row)
    stage_values = list(dict.fromkeys(str(row["stage"]) for row in sample_rows))
    stage_palette = ["#0f766e", "#2563eb", "#7c3aed", "#db2777", "#ea580c", "#65a30d", "#0891b2", "#9333ea", "#be123c", "#475569"]
    stage_colors = {stage: stage_palette[index % len(stage_palette)] for index, stage in enumerate(stage_values)}
    (figure_dir / "global_pca_by_study.svg").write_text(
        scatter_svg(global_plot_rows, "study", STUDY_COLORS, "Global exploratory PCA by study", "Top variable genes; gene-wise standardized log2(TPM + 1). Descriptive only.", tuple(global_explained)),
        encoding="utf-8",
    )
    (figure_dir / "global_pca_by_stage.svg").write_text(
        scatter_svg(global_plot_rows, "stage", stage_colors, "Global exploratory PCA by stage", "Study and biological context may be confounded; this is not a batch-effect test.", tuple(global_explained)),
        encoding="utf-8",
    )
    for study in studies:
        block = expression[:, study_column_slices[study]]
        coordinates, explained = pca_coordinates(block, int(config["pca_top_variable_genes"]))
        selected_samples = [row for row in sample_rows if row["study"] == study]
        plot_rows: list[dict[str, object]] = []
        conditions = list(dict.fromkeys(str(row["condition"]) for row in selected_samples))
        condition_colors = {condition: stage_palette[index % len(stage_palette)] for index, condition in enumerate(conditions)}
        for index, sample in enumerate(selected_samples):
            row = {
                "scope": study,
                "study": study,
                "sample_id": sample["sample_id"],
                "sample_key": sample["sample_key"],
                "condition": sample["condition"],
                "stage": sample["stage"],
                "PC1": float(coordinates[index, 0]),
                "PC2": float(coordinates[index, 1]),
                "PC1_variance_percent": float(explained[0] * 100),
                "PC2_variance_percent": float(explained[1] * 100),
            }
            pca_rows.append(row)
            plot_rows.append(row)
        (figure_dir / f"{study}_pca_by_condition.svg").write_text(
            scatter_svg(plot_rows, "condition", condition_colors, f"{study} PCA by condition", "Within-study log2(TPM + 1) PCA; descriptive view of biological samples.", tuple(explained)),
            encoding="utf-8",
        )
    pca_fields = ["scope", "study", "sample_id", "sample_key", "condition", "stage", "PC1", "PC2", "PC1_variance_percent", "PC2_variance_percent"]
    write_tsv(data_dir / "pca_coordinates.tsv", pca_fields, pca_rows)

    transformed = np.log2(expression + 1.0)
    correlation = np.corrcoef(transformed.T)
    correlation_rows = []
    for row_index, sample_key in enumerate(sample_keys):
        correlation_rows.append({"sample_key": sample_key, **{key: f"{correlation[row_index, column_index]:.6f}" for column_index, key in enumerate(sample_keys)}})
    write_tsv(data_dir / "sample_correlation.tsv", ["sample_key", *sample_keys], correlation_rows)
    correlation_labels = [f'{row["study"]}:{row["sample_id"]}' for row in sample_rows]
    (figure_dir / "sample_correlation.svg").write_text(
        heatmap_svg(correlation, correlation_labels, correlation_labels, "Sample correlation", "Pearson correlation of log2(TPM + 1); ordered by study and frozen sample order.", -1.0, 1.0),
        encoding="utf-8",
    )

    candidate_source = root / config["selected_genes"]
    curated_records: list[dict[str, str]] = []
    candidate_metadata: dict[str, dict[str, str]] = {}
    if candidate_source.suffix.lower() == ".tsv":
        candidate_genes, candidate_groups, candidate_metadata, curated_records = parse_curated_candidate_genes(candidate_source)
        shutil.copyfile(candidate_source, data_dir / "histone_modifier_gene_set_v2.tsv")
    else:
        candidate_genes, candidate_groups = parse_candidate_genes(candidate_source)
    selected_rows = []
    candidate_resolution: dict[str, str | None] = {}
    candidate_resolution_status: Counter[str] = Counter()
    for gene in candidate_genes:
        memberships = [group for group, genes in candidate_groups.items() if gene in genes]
        alternatives = sorted(set(previous_id_index.get(gene, [])))
        curated = candidate_metadata.get(gene)
        if curated:
            resolved_gene = curated["resolved_gene_id"]
            resolution_status = curated["mapping_status"]
            if resolved_gene not in gene_index:
                raise ValueError(f"curated candidate is absent from atlas gene universe: {gene} -> {resolved_gene}")
            if resolution_status == "EXACT_SM_V10" and resolved_gene != gene:
                raise ValueError(f"invalid exact curated mapping: {gene} -> {resolved_gene}")
            if resolution_status == "UNIQUE_PREVIOUS_STABLE_ID" and alternatives != [resolved_gene]:
                raise ValueError(f"curated previous-ID mapping conflicts with reference: {gene} -> {resolved_gene}")
        elif gene in gene_index:
            resolved_gene = gene
            resolution_status = "CURRENT_ID"
        elif len(alternatives) == 1:
            resolved_gene = alternatives[0]
            resolution_status = "RESOLVED_PREVIOUS_STABLE_ID"
        elif alternatives:
            resolved_gene = None
            resolution_status = "AMBIGUOUS_PREVIOUS_STABLE_ID"
        else:
            resolved_gene = None
            resolution_status = "NOT_FOUND_WBPS19"
        candidate_resolution[gene] = resolved_gene
        candidate_resolution_status[resolution_status] += 1
        annotation = annotations.get(resolved_gene or gene, {})
        selected_rows.append({
            "gene_id": gene,
            "resolved_gene_id": resolved_gene or "",
            "resolution_status": resolution_status,
            "alternative_gene_ids": ";".join(alternatives),
            "display_label": annotation.get("display_label", gene),
            "functional_name": annotation.get("functional_name", ""),
            "biotype": annotation.get("biotype", ""),
            "source_accession": annotation.get("source_accession", ""),
            "groups": ";".join(memberships),
            "present_in_reference": "YES" if resolved_gene else "NO",
            **{field: curated.get(field, "") if curated else "" for field in (
                "scope_tiers", "broad_classes", "functional_groups", "curation_statuses",
                "record_ids", "evidence_sources", "source_urls", "curation_notes",
            )},
        })
    write_tsv(
        data_dir / "selected_genes.tsv",
        ["gene_id", "resolved_gene_id", "resolution_status", "alternative_gene_ids", "display_label", "functional_name", "biotype", "source_accession", "groups", "present_in_reference", "scope_tiers", "broad_classes", "functional_groups", "curation_statuses", "record_ids", "evidence_sources", "source_urls", "curation_notes"],
        selected_rows,
    )

    # Differential-expression values remain exactly those emitted by each study.
    contrast_rows: list[dict[str, object]] = []
    contrast_keys: list[str] = []
    contrast_index: dict[str, int] = {}
    lfc_columns: list[np.ndarray] = []
    padj_columns: list[np.ndarray] = []
    alpha = float(config["significance"]["alpha"])
    lfc_threshold = float(config["significance"]["absolute_log2_fold_change"])
    study_contrast_counters: Counter[str] = Counter()
    for study_entry in config["studies"]:
        if not study_entry["differential_expression"]:
            continue
        study = study_entry["id"]
        spec_path, spec = study_spec(root, study_entry)
        results_path = root / "results" / study / "differential_expression/DEGs_all_results.tsv"
        input_files.append(spec_path)
        specs = {item["id"]: item for item in spec["contrasts"]}
        lfc_for_study = {key: np.full(len(all_genes), np.nan) for key in specs}
        padj_for_study = {key: np.full(len(all_genes), np.nan) for key in specs}
        if "contrast_sources" in study_entry:
            result_sources = approved_contrast_sources(root, study_entry, spec)
            input_files.extend(path for _, path in result_sources)
        elif results_path.is_file():
            input_files.append(results_path)
            result_sources = [(None, results_path)]
        else:
            contrast_dir = root / "results" / study / "differential_expression/contrasts"
            result_sources = [
                (contrast_id, contrast_dir / f"DEG_{contrast_id}.tsv")
                for contrast_id in specs
            ]
            missing_sources = [path for _, path in result_sources if not path.is_file()]
            if missing_sources:
                raise ValueError(
                    f"missing differential-expression source for {study}: "
                    f"{missing_sources[:3]}"
                )
            input_files.extend(path for _, path in result_sources)
        for frozen_contrast, source_path in result_sources:
            for row in read_tsv(source_path):
                contrast = frozen_contrast or row["contrast"]
                gene = row["gene_id"]
                if contrast not in specs or gene not in gene_index:
                    continue
                index = gene_index[gene]
                lfc_for_study[contrast][index] = safe_float(row.get("log2FoldChange"))
                padj_for_study[contrast][index] = safe_float(row.get("padj"))
        for contrast in spec["contrasts"]:
            contrast_id = contrast["id"]
            key = f"{study}:{contrast_id}"
            study_contrast_counters[study] += 1
            display_id = f"{study}:C{study_contrast_counters[study]:02d}"
            lfc = lfc_for_study[contrast_id]
            padj = padj_for_study[contrast_id]
            tested = np.isfinite(lfc)
            significant = tested & np.isfinite(padj) & (padj < alpha) & (np.abs(lfc) >= lfc_threshold)
            contrast_index[key] = len(contrast_keys)
            contrast_keys.append(key)
            lfc_columns.append(lfc)
            padj_columns.append(padj)
            contrast_rows.append({
                "contrast_key": key,
                "display_id": display_id,
                "study": study,
                "factor": contrast["factor"],
                "numerator": contrast["numerator"],
                "denominator": contrast["denominator"],
                "direction": contrast["direction"],
                "description": contrast["description"],
                "tested_genes": int(tested.sum()),
                "significant_genes": int(significant.sum()),
                "upregulated": int((significant & (lfc > 0)).sum()),
                "downregulated": int((significant & (lfc < 0)).sum()),
                "alpha": alpha,
                "absolute_log2_fold_change": lfc_threshold,
                "interpretation_warning": study_entry.get("interpretation_warning", ""),
                "analysis_id": spec.get("analysis_id", ""),
            })
    lfc_matrix = np.column_stack(lfc_columns)
    padj_matrix = np.column_stack(padj_columns)
    contrast_fields = ["contrast_key", "display_id", "study", "factor", "numerator", "denominator", "direction", "description", "tested_genes", "significant_genes", "upregulated", "downregulated", "alpha", "absolute_log2_fold_change", "interpretation_warning", "analysis_id"]
    write_tsv(data_dir / "contrast_catalog.tsv", contrast_fields, contrast_rows)
    # Explicit NA values keep the rectangular TSV contract unambiguous while
    # avoiding trailing empty fields that Git interprets as whitespace errors.
    lfc_rows = [{"gene_id": gene, **{key: "NA" if not math.isfinite(lfc_matrix[index, col]) else f"{lfc_matrix[index, col]:.8g}" for col, key in enumerate(contrast_keys)}} for index, gene in enumerate(all_genes)]
    padj_rows = [{"gene_id": gene, **{key: "NA" if not math.isfinite(padj_matrix[index, col]) else f"{padj_matrix[index, col]:.8g}" for col, key in enumerate(contrast_keys)}} for index, gene in enumerate(all_genes)]
    write_tsv(data_dir / "gene_contrast_log2fc.tsv", ["gene_id", *contrast_keys], lfc_rows)
    write_tsv(data_dir / "gene_contrast_padj.tsv", ["gene_id", *contrast_keys], padj_rows)

    # Per-study expression summaries preserve interpretability without fitting a cross-study model.
    expression_summary_rows: list[dict[str, object]] = []
    for study in studies:
        block = expression[:, study_column_slices[study]]
        for index, gene in enumerate(all_genes):
            values = block[index, :]
            expression_summary_rows.append({
                "gene_id": gene,
                "functional_name": annotations[gene]["functional_name"],
                "study": study,
                "n_samples": block.shape[1],
                "mean_tpm": f"{float(np.mean(values)):.8g}",
                "median_tpm": f"{float(np.median(values)):.8g}",
                "max_tpm": f"{float(np.max(values)):.8g}",
                "fraction_samples_tpm_gt_1": f"{float(np.mean(values > 1)):.8g}",
            })
    write_tsv(data_dir / "gene_expression_summary.tsv", ["gene_id", "functional_name", "study", "n_samples", "mean_tpm", "median_tpm", "max_tpm", "fraction_samples_tpm_gt_1"], expression_summary_rows)

    present_candidates = list(dict.fromkeys(
        resolved
        for gene in candidate_genes
        if (resolved := candidate_resolution[gene]) is not None
    ))
    candidate_indices = [gene_index[gene] for gene in present_candidates]
    requested_by_resolved: dict[str, list[str]] = defaultdict(list)
    for requested_gene, resolved_gene in candidate_resolution.items():
        if resolved_gene:
            requested_by_resolved[resolved_gene].append(requested_gene)
    candidate_labels = []
    for gene in present_candidates:
        functional_name = annotations[gene]["functional_name"]
        if len(functional_name) > 58:
            functional_name = functional_name[:55].rstrip() + "…"
        legacy_ids = [item for item in requested_by_resolved[gene] if item != gene]
        legacy_label = f" (legacy: {', '.join(legacy_ids)})" if legacy_ids else ""
        candidate_labels.append(f"{gene}{legacy_label} — {functional_name}")
    resolved_candidate_groups = {
        group: list(dict.fromkeys(
            resolved
            for requested in genes
            if (resolved := candidate_resolution.get(requested)) is not None
        ))
        for group, genes in candidate_groups.items()
    }
    input_files.extend(build_gene_collections(root, output, annotations, set(all_genes),
                                             resolved_candidate_groups, candidate_resolution, config["reference_id"]))
    chip_page = output / "chipseq.html"
    if chip_page.exists():
        chip_page.write_text(version_report_assets(chip_page.read_text(encoding="utf-8"), output), encoding="utf-8")
    selected_lfc = lfc_matrix[candidate_indices, :]
    selected_padj = padj_matrix[candidate_indices, :]
    significant = np.isfinite(selected_padj) & (selected_padj < alpha) & (np.abs(selected_lfc) >= lfc_threshold)
    contrast_display_ids = [str(row["display_id"]) for row in contrast_rows]
    (figure_dir / "selected_gene_contrast_matrix.svg").write_text(
        heatmap_svg(selected_lfc, candidate_labels, contrast_display_ids, "Selected-gene effect matrix", "Color is original per-study DESeq2 log2 fold change; dot marks the frozen significance rule. Full orientations are in the contrast catalog.", -5.0, 5.0, significant),
        encoding="utf-8",
    )

    # Candidate expression by study-condition context.
    contexts: list[str] = []
    context_columns: list[list[int]] = []
    context_metadata: list[dict[str, object]] = []
    for study in studies:
        study_samples = [(index, row) for index, row in enumerate(sample_rows) if row["study"] == study]
        for condition in dict.fromkeys(str(row["condition"]) for _, row in study_samples):
            matching = [(index, row) for index, row in study_samples if row["condition"] == condition]
            contexts.append(f"{study}:{condition}")
            context_columns.append([index for index, _ in matching])
            context_metadata.append(matching[0][1])
    context_matrix = np.column_stack([np.mean(np.log2(expression[candidate_indices][:, indices] + 1.0), axis=1) for indices in context_columns])
    row_mean = np.mean(context_matrix, axis=1, keepdims=True)
    row_sd = np.std(context_matrix, axis=1, keepdims=True)
    row_sd[row_sd == 0] = 1.0
    context_z = (context_matrix - row_mean) / row_sd
    (figure_dir / "selected_gene_expression_heatmap.svg").write_text(
        heatmap_svg(context_z, candidate_labels, contexts, "Selected-gene expression landscape", "Gene-wise z-score of mean log2(TPM + 1) in each study-condition context.", -2.5, 2.5),
        encoding="utf-8",
    )
    for study in studies:
        columns = [index for index, context in enumerate(contexts) if context.startswith(f"{study}:")]
        columns.sort(key=lambda index: life_cycle_sort_key(context_metadata[index]))
        study_matrix = context_matrix[:, columns]
        study_mean = np.mean(study_matrix, axis=1, keepdims=True)
        study_sd = np.std(study_matrix, axis=1, keepdims=True)
        study_sd[study_sd == 0] = 1.0
        study_z = (study_matrix - study_mean) / study_sd
        cluster = average_linkage_cluster(study_z)
        study_labels = [f"{gene} — {annotations[gene]['functional_name'][:36].rstrip()}" for gene in present_candidates]
        condition_labels = [contexts[index].split(":", 1)[1] for index in columns]
        (figure_dir / f"{study}_candidate_expression_heatmap.svg").write_text(
            heatmap_svg(
                study_z, study_labels, condition_labels,
                f"{study} epigenetic-gene expression",
                "Within-study gene-wise z-score of condition-mean log2(TPM + 1); genes clustered, conditions in life-cycle and time order where annotated.",
                -2.5, 2.5, row_cluster=cluster,
            ),
            encoding="utf-8",
        )

    contrast_labels = contrast_display_ids
    (figure_dir / "contrast_significance_summary.svg").write_text(
        bar_svg(contrast_labels, [
            ("Upregulated", [float(row["upregulated"]) for row in contrast_rows], "#dc2626"),
            ("Downregulated", [float(row["downregulated"]) for row in contrast_rows], "#2563eb"),
        ], "Differential-expression overview", "Original per-study DESeq2 contrasts; padj < 0.05 and |log2FC| ≥ 1.", "Significant genes"),
        encoding="utf-8",
    )

    (figure_dir / "study_composition.svg").write_text(
        bar_svg(studies, [
            ("Biological samples", [float(row["biological_samples"]) for row in study_rows], "#2563eb"),
            ("Technical runs", [float(row["technical_runs"]) for row in study_rows], "#94a3b8"),
        ], "Study composition", "Technical runs are aggregated before inference; biological samples are the atlas observations.", "Count"),
        encoding="utf-8",
    )

    qc_labels = studies
    (figure_dir / "qc_classification.svg").write_text(
        bar_svg(qc_labels, [
            ("PASS", [float(row["qc_pass"]) for row in study_rows], QC_COLORS["PASS"]),
            ("REVIEW", [float(row["qc_review"]) for row in study_rows], QC_COLORS["REVIEW"]),
            ("Not recorded per sample", [float(row["qc_not_recorded"]) for row in study_rows], QC_COLORS["NOT_RECORDED"]),
        ], "QC classification", "Absence of a public per-sample QC table is shown explicitly and is not interpreted as failure.", "Biological samples"),
        encoding="utf-8",
    )

    # Compact browser payload. Arrays avoid repeating field names for every value.
    study_profiles = []
    for entry in config["studies"]:
        study = entry["id"]
        spec_path, spec = study_spec(root, entry)
        if spec_path.is_file():
            input_files.append(spec_path)
        qc_path = root / entry["qc_table"]
        reasons = {row["sample_id"]: row.get("reason", "") for row in read_tsv(qc_path)}
        profile = {
            "study": study,
            "reference_id": config["reference_id"],
            "differential_expression": entry["differential_expression"],
            "design": spec.get("design", {}),
            "filter": spec.get("filter", {}),
            "parameters": spec.get("parameters", {}),
            "provider": spec.get("provider", ""),
            "test": spec.get("test", ""),
            "analysis_id": spec.get("analysis_id", ""),
            "blocked_reason": spec.get("reason", ""),
            "interpretation_warning": entry.get("interpretation_warning", ""),
            "samples": [{**row, "qc_reason": reasons.get(str(row["sample_id"]), "")} for row in sample_rows if row["study"] == study],
            "contrasts": [row for row in contrast_rows if row["study"] == study],
            "sources": [{"path": path.relative_to(root).as_posix(), "sha256": sha256(path)} for path in sorted(set(input_files)) if study in path.relative_to(root).parts and "hf_chipseq" not in path.relative_to(root).parts],
        }
        study_profiles.append(profile)
        (data_dir / f"study_{study}.json").write_text(json.dumps(profile, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    payload = {
        "study_profiles": study_profiles,
        "schema": "1.0",
        "reference": config["reference_id"],
        "samples": [[row[field] for field in (
            "study", "sample_id", "condition", "stage", "tissue", "sex",
            "infection_mode", "treatment", "time_hours", "batch", "replicate",
            "qc_classification",
        )] for row in sample_rows],
        "sample_fields": [
            "study", "sample_id", "condition", "stage", "tissue", "sex",
            "infection_mode", "treatment", "time_hours", "batch", "replicate",
            "qc_classification",
        ],
        "contrasts": [[row[field] for field in ("contrast_key", "study", "numerator", "denominator", "interpretation_warning", "analysis_id")] for row in contrast_rows],
        "contrast_fields": ["contrast_key", "study", "numerator", "denominator", "interpretation_warning", "analysis_id"],
        "genes": all_genes,
        "gene_annotation_fields": ["display_label", "functional_name", "annotation_name", "description", "biotype", "description_source", "source_accession", "previous_stable_id", "chromosome", "start", "end", "strand"],
        "gene_annotations": [[annotations[gene][field] for field in ("display_label", "functional_name", "annotation_name", "description", "biotype", "description_source", "source_accession", "previous_stable_id", "chromosome", "start", "end", "strand")] for gene in all_genes],
        "gene_model_fields": ["transcript_id", "chromosome", "start", "end", "strand", "exons"],
        "gene_models": [[[model[field] for field in ("transcript_id", "chromosome", "start", "end", "strand", "exons")] for model in gene_models[gene]] for gene in all_genes],
        "tpm": [[compact_number(float(value), 5) for value in expression[index, :]] for index in range(len(all_genes))],
        "lfc": [[compact_number(float(value), 17) for value in lfc_matrix[index, :]] for index in range(len(all_genes))],
        "padj": [[compact_number(float(value), 17) for value in padj_matrix[index, :]] for index in range(len(all_genes))],
        "selected_groups": resolved_candidate_groups,
        "selected_gene_fields": [
            "gene_id", "resolved_gene_id", "resolution_status",
            "alternative_gene_ids", "functional_name", "groups", "scope_tiers",
            "broad_classes", "functional_groups", "curation_statuses", "record_ids",
            "evidence_sources", "source_urls", "curation_notes",
        ],
        "selected_gene_status": [[row[field] for field in (
            "gene_id", "resolved_gene_id", "resolution_status",
            "alternative_gene_ids", "functional_name", "groups", "scope_tiers",
            "broad_classes", "functional_groups", "curation_statuses", "record_ids",
            "evidence_sources", "source_urls", "curation_notes",
        )] for row in selected_rows],
        "external_annotation_fields": ["gene_index", "evidence_type", "stage_or_context", "assertion", "source_location", "source_record_id", "original_id", "mapping_method", "notes"],
        "external_annotation": [[gene_index[row["current_gene_id"]], *[row[field] for field in ("evidence_type", "stage_or_context", "assertion", "source_location", "source_record_id", "original_id", "mapping_method", "notes")]] for row in external_annotation_rows if row["current_gene_id"]],
        "cell_marker_fields": ["gene_index", "cell_population", "lineage_or_state", "life_stage", "marker_class", "reported_statistic", "assertion", "source_location", "source_record_id", "original_id", "mapping_method", "notes"],
        "cell_marker_source_count": len(cell_marker_rows),
        "cell_marker_mapped_gene_count": len({row["current_gene_id"] for row in cell_marker_rows if row["current_gene_id"]}),
        "cell_markers": [[gene_index[row["current_gene_id"]], *[row[field] for field in ("cell_population", "lineage_or_state", "life_stage", "marker_class", "reported_statistic", "assertion", "source_location", "source_record_id", "original_id", "mapping_method", "notes")]] for row in cell_marker_rows if row["current_gene_id"]],
        "lncrna_source_count": len(lncrna_rows),
        "lncrna_gene_link_fields": ["gene_index", "original_id", "lncrna_class", "robust_candidate", "stage_differential_expression", "associated_original_gene_id", "association_mapping_method", "source_record_id", "notes"],
        "lncrna_gene_links": [[gene_index[row["current_gene_id"]], *[row[field] for field in ("original_id", "lncrna_class", "robust_candidate", "stage_differential_expression", "associated_original_gene_id", "association_mapping_method", "source_record_id", "notes")]] for row in lncrna_gene_links],
        "almeida_source_count": len(legacy_transcript_rows),
        "almeida_evidence_fields": ["gene_index", "evidence_type", "original_id", "est_count", "published_annotation", "assertion", "source_location", "source_record_id", "mapping_method", "notes"],
        "almeida_evidence": [[gene_index[row["current_gene_id"]], *[row[field] for field in ("evidence_type", "original_id", "est_count", "published_annotation", "assertion", "source_location", "source_record_id", "mapping_method", "notes")]] for row in almeida_mapped_rows],
        "trajectory_group_by": {
            entry["id"]: entry.get("trajectory_group_by", "auto")
            for entry in config["studies"]
        },
        "chip_integration_status": chip_config["status"] if chip_config else "NOT_AVAILABLE",
        "locus_display": config["locus_display"],
        "chip_study_profiles": [{key: value for key, value in profile.items() if key not in {"source_files", "samples", "publication", "native_execution"}} for profile in chip_profiles],
        "chip_scope": "all_atlas_genes",
        "chip_mark_stage_fields": METADATA_FIELDS,
        "chip_mark_stages": [[row[field] for field in METADATA_FIELDS] for row in chip_metadata_rows],
        "chip_evidence_fields": ["gene_index", *[field for field in EVIDENCE_FIELDS if field not in {"gene_id", "functional_name"}]],
        "chip_evidence": [[gene_index[str(row["gene_id"])], *[row[field] for field in EVIDENCE_FIELDS if field not in {"gene_id", "functional_name"}]] for row in chip_evidence_rows if row["study_id"] not in delayed_studies],
        "chip_locus_fields": ["gene_index", *[field for field in LOCUS_FIELDS if field != "gene_id"]],
        "chip_loci": [[gene_index[str(row["gene_id"])], *[row[field] for field in LOCUS_FIELDS if field != "gene_id"]] for row in chip_locus_rows if row["study_id"] not in delayed_studies],
        "chip_summaries": [r for r in chip_data["lazy"]["gene_summaries"] if chip_metadata_rows[r[1]]["study_id"] in delayed_studies],
        "chip_summary_fields": chip_data["lazy"]["summary_fields"],
        "chip_category_values": chip_data["lazy"]["category_values"],
        "chip_gene_store": chip_data["lazy"]["gene_store"],
        "stage_compatibility_fields": ["study", "rna_stage", "chip_stage", "match_status"],
        "stage_compatibility": [[row[field] for field in ("study", "rna_stage", "chip_stage", "match_status")] for row in stage_compatibility_rows],
        "significance": {"alpha": alpha, "absolute_log2_fold_change": lfc_threshold},
    }
    (data_dir / "atlas_payload.js").write_text("window.HF_ATLAS=" + json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + ";\n", encoding="utf-8")

    template_path = root / "analyses/atlas/rnaseq_atlas_template.html"
    input_files.append(template_path)
    template = template_path.read_text(encoding="utf-8")
    study_table = table_html(study_fields, study_rows)
    contrast_table = table_html(["display_id", "study", "numerator", "denominator", "tested_genes", "significant_genes", "upregulated", "downregulated", "interpretation_warning"], contrast_rows)
    chip_inventory_table = table_html(["study_id", "stage_or_condition", "mark_or_factor", "n_samples", "n_controls", "qc_status"], chip_metadata_rows)
    matched_stage_rows = [row for row in stage_compatibility_rows if row["match_status"] != "NOT_COVERED"]
    stage_compatibility_table = table_html(["study", "rna_stage", "chip_stage", "match_status"], matched_stage_rows)
    chip_gene_count = len({str(row["gene_id"]) for row in chip_evidence_rows})
    chip_sample_count = sum(int(row["n_samples"]) for row in chip_metadata_rows)
    chip_stage_count = len({str(row["stage_or_condition"]) for row in chip_metadata_rows})
    chip_mark_count = len({str(row["mark_or_factor"]) for row in chip_metadata_rows})
    covered_study_count = len({str(row["study"]) for row in matched_stage_rows})
    lncrna_class_counts = Counter(str(row["lncrna_class"]) for row in lncrna_rows)
    lncrna_robust_count = sum(row["robust_candidate"] == "YES" for row in lncrna_rows)
    lncrna_de_count = sum(row["stage_differential_expression"] == "YES" for row in lncrna_rows)
    lncrna_linked_gene_count = len({str(row["current_gene_id"]) for row in lncrna_gene_links})
    almeida_type_counts = Counter(str(row["evidence_type"]) for row in legacy_transcript_rows)
    almeida_mapped_gene_count = len({str(row["current_gene_id"]) for row in almeida_mapped_rows})
    report = build_html(template, {
        "ATLAS_TITLE": html.escape(config["title"]),
        "STUDY_COUNT": str(len(studies)),
        "SAMPLE_COUNT": str(len(sample_rows)),
        "RUN_COUNT": str(sum(int(row["technical_runs"]) for row in sample_rows)),
        "GENE_COUNT": str(len(all_genes)),
        "SELECTED_GENE_COUNT": str(len(present_candidates)),
        "CONTRAST_COUNT": str(len(contrast_rows)),
        "STUDY_TABLE": study_table,
        "CONTRAST_TABLE": contrast_table,
        "CHIP_STUDY_COUNT": str(len(chip_profiles)),
        "CHIP_SEGMENT_COUNT": str(len(chip_all_segments)),
        "CHIP_STUDY_LINKS": " ".join(f'<a href="chipseq_{p["study_id"]}.html">{p["study_id"]}: design, QC &amp; provenance →</a>' for p in chip_profiles),
        "CHIP_SAMPLE_COUNT": str(chip_sample_count),
        "CHIP_STAGE_COUNT": str(chip_stage_count),
        "CHIP_MARK_COUNT": str(chip_mark_count),
        "CHIP_GENE_COUNT": str(chip_gene_count),
        "CHIP_EVIDENCE_COUNT": str(len(chip_evidence_rows)),
        "CHIP_COVERED_STUDIES": str(covered_study_count),
        "CHIP_STATUS": html.escape(chip_config["status"] if chip_config else "NOT_AVAILABLE"),
        "CHIP_INVENTORY_TABLE": chip_inventory_table,
        "STAGE_COMPATIBILITY_TABLE": stage_compatibility_table,
        "VASCONCELOS_LNCRNA_COUNT": f"{len(lncrna_rows):,}",
        "VASCONCELOS_LINC_COUNT": f"{lncrna_class_counts['LINC_RNA']:,}",
        "VASCONCELOS_ANTISENSE_COUNT": f"{lncrna_class_counts['ANTISENSE_LNCRNA']:,}",
        "VASCONCELOS_ROBUST_COUNT": f"{lncrna_robust_count:,}",
        "VASCONCELOS_DE_COUNT": f"{lncrna_de_count:,}",
        "VASCONCELOS_LINKED_GENE_COUNT": f"{lncrna_linked_gene_count:,}",
        "ALMEIDA_RECORD_COUNT": f"{len(legacy_transcript_rows):,}",
        "ALMEIDA_CONFIRMED_COUNT": f"{almeida_type_counts['PUBLISHED_ADULT_MALE_TRANSCRIPTION_SUPPORT']:,}",
        "ALMEIDA_FIRST_COUNT": f"{almeida_type_counts['PUBLISHED_FIRST_ADULT_MALE_TRANSCRIPTION_SUPPORT']:,}",
        "ALMEIDA_CONTIG_COUNT": f"{almeida_type_counts['PUBLISHED_NOVEL_TRANSCRIBED_CONTIG']:,}",
        "ALMEIDA_MAPPED_GENE_COUNT": f"{almeida_mapped_gene_count:,}",
        "REFERENCE_ID": html.escape(config["reference_id"]),
        "BATCH_STATUS": html.escape(config["batch_effect_assessment"]),
    })
    report = version_report_assets(report, output)
    atlas_report = report.replace('<body>', '<body class="atlas-page" data-portal="atlas">', 1)
    (output / "atlas.html").write_text(atlas_report, encoding="utf-8")
    (output / "machinery.html").write_text('<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Gene collection · HelixForge</title></head><body><p>The machinery is now a shared gene collection. <a href="atlas.html?collection=machinery#collection">Open the collection</a>.</p><script>const p=new URLSearchParams(location.search);p.set(\'collection\',\'machinery\');location.replace(\'atlas.html?\'+p+(location.hash&&location.hash!==\'#groups\'?location.hash:\'#collection\'));</script></body></html>\n', encoding="utf-8")
    for study in studies:
        study_report = report.replace('<body>', f'<body class="study-page" data-study="{study}">')
        nav_start, nav_end = study_report.index('  <nav '), study_report.index('  </nav>') + len('  </nav>')
        study_report = study_report[:nav_start] + '<nav aria-label="Study navigation"><a href="atlas.html#overview">← All studies</a><a href="atlas.html">Atlas home</a><a href="atlas.html#genes">Gene explorer</a></nav>' + study_report[nav_end:]
        study_report = study_report.replace('<title>', f'<title>{study} · ')
        (output / f"study_{study}.html").write_text(study_report, encoding="utf-8")

    (output / "README.md").write_text(
        "# HelixForge-SMansoni RNA-seq Atlas\n\n"
        f"Open `atlas.html` for RNA-seq or `chipseq.html` for ChIP-seq. Both explorers share gene collections, including the {len(present_candidates)}-gene machinery preset. `machinery.html` redirects to that collection. The atlas integrates {len(studies)} accepted studies, "
        f"{len(sample_rows)} biological samples, {sum(int(row['technical_runs']) for row in sample_rows)} technical runs, "
        f"{len(all_genes):,} shared genes, and {len(contrast_rows)} original per-study DESeq2 contrasts.\n\n"
        f"The reviewed HF ChIP-seq layer includes {len(chip_profiles)} study, {chip_sample_count} deposited IP libraries, "
        f"{len(chip_all_segments)} regional candidates and {chip_gene_count} genes with native associations. "
        "Gene associations and loci are available throughout the atlas gene universe. HOMER nearest-promoter assignments include source-classified intergenic regions; location class and nearest-gene assignment remain distinct.\n\n"
        "The gene explorer includes an SM_V10 coordinate-aware locus view of every native gene–segment association. "
        "The display is descriptive and does not imply sample pairing or causality.\n\n"
        f"Separate literature layers contain {len(external_annotation_rows):,} Anderson annotation/source claims and "
        f"{len(cell_marker_rows)} Wang class-dependent sporocyst stem-cell transcripts. "
        "These are published evidence, not extra bulk samples, atlas TPM values or DESeq2 contrasts. "
        "All source records, including unmapped IDs, are available in the external evidence TSVs.\n\n"
        f"The literature annotation layer additionally contains {len(lncrna_rows):,} Vasconcelos lncRNA models "
        f"({lncrna_robust_count:,} robust candidates; {lncrna_de_count:,} with published stage-DE evidence) and "
        f"{len(legacy_transcript_rows):,} Almeida historical transcript claims. Vasconcelos coordinates remain "
        "labelled as v5.2 source coordinates, and Almeida contigs without a stable gene ID remain unmapped.\n\n"
        "The gene explorer shows the current SM_V10/WBPS19 annotation for every atlas gene. "
        "Historical model claims cannot be confirmed from a stable-ID match alone.\n\n"
        "The atlas is descriptive. It does not fit a joint differential-expression model, "
        "perform batch correction, or reinterpret the original contrast orientations. "
        "Machine-readable tables are under `data/`, exportable figures under `figures/`, "
        "and all generated artifacts are recorded in `manifest.json`.\n",
        encoding="utf-8",
    )

    from chipseq_lazy import compress_downloads
    input_files.extend(apply_atlas_design(root, output))
    if delayed_studies:compress_downloads(output)
    chip_payload_bytes = len(json.dumps({key: value for key, value in payload.items() if key.startswith("chip_")}, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))
    size_report = {"scope": "static_bundle_with_gzip_downloads" if delayed_studies else "uncompressed_static_bundle", "payload_bytes": (data_dir / "atlas_payload.js").stat().st_size,
                   "chip_payload_bytes": chip_payload_bytes, "chip_studies": len(chip_profiles),
                   "chip_associations": len(chip_locus_rows), "all_gene_scope": True,
                   "chip_explorer_json_bytes": (data_dir / "chipseq/explorer.json").stat().st_size,
                   "chip_explorer_library_peaks": len(chip_data["library_peaks"]),
                   "chip_explorer_initial_asset_bytes": sum((output / name).stat().st_size for name in ["chipseq.html", "chipseq.css", "chipseq_logic.js", "chipseq_tracks.js", "chipseq.js", "data/chipseq/explorer.json", "gene_collections.js", "gene_collections.css", "data/gene_collections.json", "chip_data.js", "atlas_design.css", "atlas_theme.js"]) + sum(path.stat().st_size for path in (output / "fonts").glob("*.woff2")),
                   "shared_collection_asset_bytes": sum((output / name).stat().st_size for name in ["gene_collections.js", "gene_collections.css", "data/gene_collections.json", "chip_data.js", "atlas_design.css", "atlas_theme.js"]) + sum(path.stat().st_size for path in (output / "fonts").glob("*.woff2")),
                   "collection_reference_genes": len(annotations),
                   "collection_preset_genes": len(present_candidates),
                   "bundle_content_bytes": sum(path.stat().st_size for path in output.rglob("*") if path.is_file() and path.name not in {"manifest.json", "size_report.json"}),
                   "note": "Sparse observed associations only; no dense all-genes × marks × stages matrix. Browser timing is measured separately."}
    (output / "size_report.json").write_text(json.dumps(size_report, indent=2) + "\n", encoding="utf-8")
    output_files = sorted(path for path in output.rglob("*") if path.is_file() and path.name != "manifest.json")
    manifest = {
        "schema_version": "1.0",
        "atlas_id": config["atlas_id"],
        "status": "PASS",
        "scope": config["statistical_scope"],
        "reference_id": config["reference_id"],
        "studies": studies,
        "biological_samples": len(sample_rows),
        "technical_runs": sum(int(row["technical_runs"]) for row in sample_rows),
        "genes": len(all_genes),
        "annotated_genes": len(annotation_rows),
        "genes_with_transcript_models": sum(bool(models) for models in gene_models.values()),
        "genes_without_transcript_models": sum(not models for models in gene_models.values()),
        "transcript_models": sum(map(len, gene_models.values())),
        "transcript_model_source_sha256": model_source["source_sha256"],
        "locus_display": config["locus_display"],
        "contrasts": len(contrast_rows),
        "selected_group_memberships": sum(len(genes) for genes in candidate_groups.values()),
        "selected_unique_genes": len(candidate_genes),
        "selected_unique_genes_present": len(present_candidates),
        "selected_source_records": len(curated_records),
        "selected_excluded_review_records": sum(row["scope_tier"] == "EXCLUDED_REVIEW" for row in curated_records),
        "selected_mapping_statuses": dict(candidate_resolution_status),
        "external_annotation_records": len(external_annotation_rows),
        "external_annotation_mapped": sum(bool(row["current_gene_id"]) for row in external_annotation_rows),
        "cell_marker_records": len(cell_marker_rows),
        "cell_marker_mapped": sum(bool(row["current_gene_id"]) for row in cell_marker_rows),
        "vasconcelos_lncrna_records": len(lncrna_rows),
        "vasconcelos_linc_records": lncrna_class_counts["LINC_RNA"],
        "vasconcelos_antisense_records": lncrna_class_counts["ANTISENSE_LNCRNA"],
        "vasconcelos_robust_candidates": lncrna_robust_count,
        "vasconcelos_stage_de_records": lncrna_de_count,
        "vasconcelos_current_gene_links": len(lncrna_gene_links),
        "vasconcelos_current_genes_linked": lncrna_linked_gene_count,
        "almeida_evidence_records": len(legacy_transcript_rows),
        "almeida_mapped_records": len(almeida_mapped_rows),
        "almeida_current_genes_mapped": almeida_mapped_gene_count,
        "almeida_confirmed_prediction_records": almeida_type_counts["PUBLISHED_ADULT_MALE_TRANSCRIPTION_SUPPORT"],
        "almeida_first_support_records": almeida_type_counts["PUBLISHED_FIRST_ADULT_MALE_TRANSCRIPTION_SUPPORT"],
        "almeida_novel_contig_records": almeida_type_counts["PUBLISHED_NOVEL_TRANSCRIBED_CONTIG"],
        "selected_previous_ids_resolved": candidate_resolution_status["RESOLVED_PREVIOUS_STABLE_ID"] + candidate_resolution_status["UNIQUE_PREVIOUS_STABLE_ID"],
        "selected_ambiguous_previous_ids": candidate_resolution_status["AMBIGUOUS_PREVIOUS_STABLE_ID"],
        "selected_unresolved_ids": candidate_resolution_status["NOT_FOUND_WBPS19"],
        "batch_effect_assessment": config["batch_effect_assessment"],
        "chipseq_integration_status": chip_config["status"] if chip_config else "NOT_AVAILABLE",
        "chipseq_studies": [profile["study_id"] for profile in chip_profiles],
        "chipseq_scope": "all_atlas_genes",
        "chipseq_controls": sum(int(row["n_controls"]) for row in chip_metadata_rows),
        "chipseq_consensus_segments": len(chip_all_segments),
        "chipseq_shared_core_segments": len(chip_data["core_segments"]),
        "chipseq_single_library_regions": sum(r["support"] == 1 for r in chip_all_segments),
        "chipseq_intergenic_segments": sum(row["category"] == "intergenic" for row in chip_all_segments),
        "chipseq_gene_peak_loci": len(chip_locus_rows),
        "chipseq_samples": chip_sample_count,
        "chipseq_stages": chip_stage_count,
        "chipseq_marks_or_factors": chip_mark_count,
        "chipseq_genes_with_evidence": chip_gene_count,
        "chipseq_gene_mark_stage_associations": len(chip_evidence_rows),
        "chipseq_gene_position_records": len(chip_position_rows),
        "rna_studies_with_stage_context": covered_study_count,
        "global_pca_interpretation": "EXPLORATORY_NOT_A_BATCH_TEST",
        "inputs": [{"path": path.relative_to(root).as_posix(), "sha256": sha256(path)} for path in sorted(set(input_files))],
        "outputs": [{"path": path.relative_to(output).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in output_files],
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({key: manifest[key] for key in ("status", "studies", "biological_samples", "technical_runs", "genes", "contrasts", "selected_unique_genes_present")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
