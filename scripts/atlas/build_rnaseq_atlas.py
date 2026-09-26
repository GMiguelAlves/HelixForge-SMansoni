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
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np


STUDY_COLORS = {
    "PRJNA602528": "#0f766e",
    "PRJNA597909": "#2563eb",
    "PRJEB14695": "#7c3aed",
    "PRJEB32839": "#db2777",
}
QC_COLORS = {"PASS": "#16a34a", "REVIEW": "#d97706", "FAIL": "#dc2626", "NOT_RECORDED": "#64748b"}


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


def compact_number(value: float, digits: int = 5) -> float | None:
    if not math.isfinite(value):
        return None
    return float(f"{value:.{digits}g}")


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
) -> str:
    cell_w = max(14, min(34, 850 // max(1, len(column_labels))))
    cell_h = max(11, min(24, 950 // max(1, len(row_labels))))
    left, top, right, bottom = 220, 125, 80, 250
    width = left + len(column_labels) * cell_w + right
    height = top + len(row_labels) * cell_h + bottom
    elements = [
        f'<text class="title" x="24" y="36">{xml_text(title)}</text>',
        f'<text class="subtitle" x="24" y="58">{xml_text(subtitle)}</text>',
    ]
    for row_index, label in enumerate(row_labels):
        y = top + row_index * cell_h
        elements.append(f'<text class="small" text-anchor="end" x="{left - 8}" y="{y + cell_h * 0.78:.1f}">{xml_text(label)}</text>')
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
        elements.append(f'<circle cx="{legend_x + 250}" cy="{legend_y + 6}" r="3" fill="#111827"/><text class="small" x="{legend_x + 260}" y="{legend_y + 10}">significant (padj &lt; 0.05 and |log2FC| ≥ 1)</text>')
    return svg_shell(width, height, "".join(elements), title, subtitle)


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
    input_files = [config_path, root / config["selected_genes"], Path(__file__).resolve()]
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
                "source_sample_name": meta.get("source_sample_name", ""),
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
                "technical_runs": int(meta.get("technical_runs", "1")),
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
    expression = np.concatenate(expression_blocks, axis=1)
    gene_index = {gene: index for index, gene in enumerate(all_genes)}
    sample_fields = [
        "study", "sample_id", "sample_key", "source_sample_name", "biosample", "condition", "stage", "tissue",
        "sex", "infection_mode", "treatment", "time_hours", "batch", "replicate", "technical_runs",
        "qc_classification", "trim_retention_percent", "salmon_mapping_percent", "reference_id",
    ]
    write_tsv(data_dir / "sample_inventory.tsv", sample_fields, sample_rows)

    study_rows: list[dict[str, object]] = []
    for study_entry in config["studies"]:
        study = study_entry["id"]
        selected = [row for row in sample_rows if row["study"] == study]
        qc_counts = Counter(str(row["qc_classification"]) for row in selected)
        mapping_values = [float(row["salmon_mapping_percent"]) for row in selected if row["salmon_mapping_percent"] is not None]
        retention_values = [float(row["trim_retention_percent"]) for row in selected if row["trim_retention_percent"] is not None]
        aggregate_qc_path = root / "results" / study / "qc/qc_summary.json"
        aggregate_qc = {}
        if aggregate_qc_path.is_file():
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
        })
    study_fields = [
        "study", "biological_samples", "technical_runs", "conditions", "stages", "qc_pass", "qc_review", "qc_fail",
        "qc_not_recorded", "mapping_mean_percent", "mapping_min_percent", "mapping_max_percent",
        "trim_retention_mean_percent", "trim_retention_min_percent", "trim_retention_max_percent", "differential_expression",
    ]
    write_tsv(data_dir / "study_summary.tsv", study_fields, study_rows)

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

    candidate_genes, candidate_groups = parse_candidate_genes(root / config["selected_genes"])
    selected_rows = []
    for gene in candidate_genes:
        memberships = [group for group, genes in candidate_groups.items() if gene in genes]
        selected_rows.append({"gene_id": gene, "groups": ";".join(memberships), "present_in_reference": "YES" if gene in gene_index else "NO"})
    write_tsv(data_dir / "selected_genes.tsv", ["gene_id", "groups", "present_in_reference"], selected_rows)

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
        spec_path = root / "config" / study / "de_spec.json"
        results_path = root / "results" / study / "differential_expression/DEGs_all_results.tsv"
        input_files.extend([spec_path, results_path])
        spec = json.loads(spec_path.read_text(encoding="utf-8"))
        specs = {item["id"]: item for item in spec["contrasts"]}
        lfc_for_study = {key: np.full(len(all_genes), np.nan) for key in specs}
        padj_for_study = {key: np.full(len(all_genes), np.nan) for key in specs}
        for row in read_tsv(results_path):
            contrast = row["contrast"]
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
            })
    lfc_matrix = np.column_stack(lfc_columns)
    padj_matrix = np.column_stack(padj_columns)
    contrast_fields = ["contrast_key", "display_id", "study", "factor", "numerator", "denominator", "direction", "description", "tested_genes", "significant_genes", "upregulated", "downregulated", "alpha", "absolute_log2_fold_change"]
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
                "study": study,
                "n_samples": block.shape[1],
                "mean_tpm": f"{float(np.mean(values)):.8g}",
                "median_tpm": f"{float(np.median(values)):.8g}",
                "max_tpm": f"{float(np.max(values)):.8g}",
                "fraction_samples_tpm_gt_1": f"{float(np.mean(values > 1)):.8g}",
            })
    write_tsv(data_dir / "gene_expression_summary.tsv", ["gene_id", "study", "n_samples", "mean_tpm", "median_tpm", "max_tpm", "fraction_samples_tpm_gt_1"], expression_summary_rows)

    present_candidates = [gene for gene in candidate_genes if gene in gene_index]
    candidate_indices = [gene_index[gene] for gene in present_candidates]
    selected_lfc = lfc_matrix[candidate_indices, :]
    selected_padj = padj_matrix[candidate_indices, :]
    significant = np.isfinite(selected_padj) & (selected_padj < alpha) & (np.abs(selected_lfc) >= lfc_threshold)
    contrast_display_ids = [str(row["display_id"]) for row in contrast_rows]
    (figure_dir / "selected_gene_contrast_matrix.svg").write_text(
        heatmap_svg(selected_lfc, present_candidates, contrast_display_ids, "Selected-gene effect matrix", "Color is original per-study DESeq2 log2 fold change; dot marks the frozen significance rule. Full orientations are in the contrast catalog.", -5.0, 5.0, significant),
        encoding="utf-8",
    )

    # Candidate expression by study-condition context.
    contexts: list[str] = []
    context_columns: list[list[int]] = []
    for study in studies:
        study_samples = [(index, row) for index, row in enumerate(sample_rows) if row["study"] == study]
        for condition in dict.fromkeys(str(row["condition"]) for _, row in study_samples):
            contexts.append(f"{study}:{condition}")
            context_columns.append([index for index, row in study_samples if row["condition"] == condition])
    context_matrix = np.column_stack([np.mean(np.log2(expression[candidate_indices][:, indices] + 1.0), axis=1) for indices in context_columns])
    row_mean = np.mean(context_matrix, axis=1, keepdims=True)
    row_sd = np.std(context_matrix, axis=1, keepdims=True)
    row_sd[row_sd == 0] = 1.0
    context_z = (context_matrix - row_mean) / row_sd
    (figure_dir / "selected_gene_expression_heatmap.svg").write_text(
        heatmap_svg(context_z, present_candidates, contexts, "Selected-gene expression landscape", "Gene-wise z-score of mean log2(TPM + 1) in each study-condition context.", -2.5, 2.5),
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
    payload = {
        "schema": "1.0",
        "reference": config["reference_id"],
        "samples": [[row[field] for field in ("study", "sample_id", "condition", "stage", "tissue", "sex", "treatment", "batch", "qc_classification")] for row in sample_rows],
        "sample_fields": ["study", "sample_id", "condition", "stage", "tissue", "sex", "treatment", "batch", "qc_classification"],
        "contrasts": [[row[field] for field in ("contrast_key", "study", "numerator", "denominator")] for row in contrast_rows],
        "contrast_fields": ["contrast_key", "study", "numerator", "denominator"],
        "genes": all_genes,
        "tpm": [[compact_number(float(value), 5) for value in expression[index, :]] for index in range(len(all_genes))],
        "lfc": [[compact_number(float(value), 5) for value in lfc_matrix[index, :]] for index in range(len(all_genes))],
        "padj": [[compact_number(float(value), 4) for value in padj_matrix[index, :]] for index in range(len(all_genes))],
        "selected_groups": candidate_groups,
        "significance": {"alpha": alpha, "absolute_log2_fold_change": lfc_threshold},
    }
    (data_dir / "atlas_payload.js").write_text("window.HF_ATLAS=" + json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + ";\n", encoding="utf-8")

    template_path = root / "analyses/atlas/rnaseq_atlas_template.html"
    input_files.append(template_path)
    template = template_path.read_text(encoding="utf-8")
    study_table = table_html(study_fields, study_rows)
    contrast_table = table_html(["display_id", "study", "numerator", "denominator", "tested_genes", "significant_genes", "upregulated", "downregulated"], contrast_rows)
    report = build_html(template, {
        "ATLAS_TITLE": html.escape(config["title"]),
        "STUDY_COUNT": str(len(studies)),
        "SAMPLE_COUNT": str(len(sample_rows)),
        "RUN_COUNT": str(sum(int(row["technical_runs"]) for row in sample_rows)),
        "GENE_COUNT": str(len(all_genes)),
        "CONTRAST_COUNT": str(len(contrast_rows)),
        "STUDY_TABLE": study_table,
        "CONTRAST_TABLE": contrast_table,
        "REFERENCE_ID": html.escape(config["reference_id"]),
        "BATCH_STATUS": html.escape(config["batch_effect_assessment"]),
    })
    (output / "atlas.html").write_text(report, encoding="utf-8")
    (output / "README.md").write_text(
        "# HelixForge-SMansoni RNA-seq Atlas\n\n"
        "Open `atlas.html` in a modern browser. The report integrates four accepted studies, "
        "128 biological samples, 318 technical runs, 9,914 shared genes, and 22 original "
        "per-study DESeq2 contrasts.\n\n"
        "The atlas is descriptive. It does not fit a joint differential-expression model, "
        "perform batch correction, or reinterpret the original contrast orientations. "
        "Machine-readable tables are under `data/`, exportable figures under `figures/`, "
        "and all generated artifacts are recorded in `manifest.json`.\n",
        encoding="utf-8",
    )

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
        "contrasts": len(contrast_rows),
        "selected_group_memberships": sum(len(genes) for genes in candidate_groups.values()),
        "selected_unique_genes": len(candidate_genes),
        "selected_unique_genes_present": len(present_candidates),
        "batch_effect_assessment": config["batch_effect_assessment"],
        "global_pca_interpretation": "EXPLORATORY_NOT_A_BATCH_TEST",
        "inputs": [{"path": path.relative_to(root).as_posix(), "sha256": sha256(path)} for path in sorted(set(input_files))],
        "outputs": [{"path": path.relative_to(output).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in output_files],
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({key: manifest[key] for key in ("status", "studies", "biological_samples", "technical_runs", "genes", "contrasts", "selected_unique_genes_present")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
