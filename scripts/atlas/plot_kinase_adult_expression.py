#!/usr/bin/env python3
"""Create focused kinase-expression figures from frozen atlas TPM matrices.

The figures are descriptive. Replicates are summarized only within an exact
biological condition; distinct stages, sexes and pairing states are never
pooled. Cross-study values are displayed together without fitting a joint
statistical model.
"""

from __future__ import annotations

import argparse
import html
import math
from pathlib import Path

import numpy as np
import pandas as pd


GENES = [
    "Smp_073340",
    "Smp_306370",
    "Smp_096640",
    "Smp_139800",
    "Smp_333950",
    "Smp_158560",
]
TARGET = "Smp_073340"

PRJEB_ORDER = [
    "egg_mixed",
    "miracidium_mixed",
    "sporocyst_1d_mixed",
    "sporocyst_5d_mixed",
    "sporocyst_32d_female",
    "sporocyst_32d_male",
    "sporocyst_32d_mixed",
    "cercaria_female",
    "cercaria_male",
    "cercaria_mixed",
    "schistosomulum_2d_female",
    "schistosomulum_2d_male",
    "schistosomulum_2d_mixed",
    "juvenile_26d_female",
    "juvenile_26d_male",
]

LABELS = {
    "egg_mixed": "Egg\nmixed",
    "miracidium_mixed": "Miracidium\nmixed",
    "sporocyst_1d_mixed": "Sporocyst 1d\nmixed",
    "sporocyst_5d_mixed": "Sporocyst 5d\nmixed",
    "sporocyst_32d_female": "Sporocyst 32d\nfemale",
    "sporocyst_32d_male": "Sporocyst 32d\nmale",
    "sporocyst_32d_mixed": "Sporocyst 32d\nmixed",
    "cercaria_female": "Cercaria\nfemale",
    "cercaria_male": "Cercaria\nmale",
    "cercaria_mixed": "Cercaria\nmixed",
    "schistosomulum_2d_female": "Schistosomulum 2d\nfemale",
    "schistosomulum_2d_male": "Schistosomulum 2d\nmale",
    "schistosomulum_2d_mixed": "Schistosomulum 2d\nmixed",
    "juvenile_26d_female": "Juvenile 26d\nfemale",
    "juvenile_26d_male": "Juvenile 26d\nmale",
}

SERIES = {
    "female_paired": ("Female · paired", "#2563eb", "circle"),
    "female_unpaired": ("Female · unpaired", "#7c3aed", "square"),
    "male_paired": ("Male · paired", "#059669", "circle"),
    "male_unpaired": ("Male · unpaired", "#d97706", "square"),
}

GONAD_CONDITIONS = [
    "ovary_mixed_sex",
    "ovary_single_sex",
    "testis_mixed_sex",
    "testis_single_sex",
]

GONAD_LABELS = {
    "ovary_mixed_sex": "Ovary\nmixed-sex",
    "ovary_single_sex": "Ovary\nsingle-sex",
    "testis_mixed_sex": "Testis\nmixed-sex",
    "testis_single_sex": "Testis\nsingle-sex",
}

GONAD_COLORS = {
    "ovary_mixed_sex": "#db2777",
    "ovary_single_sex": "#f472b6",
    "testis_mixed_sex": "#2563eb",
    "testis_single_sex": "#60a5fa",
}


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def load_study(root: Path, study: str) -> pd.DataFrame:
    tpm_path = root / "results" / study / "expression" / "tpm_matrix.tsv"
    metadata_path = root / "metadata" / study / "samples.tsv"
    tpm = pd.read_csv(tpm_path, sep="\t", index_col=0)
    metadata = pd.read_csv(metadata_path, sep="\t").set_index("sample_id")
    missing = [gene for gene in GENES if gene not in tpm.index]
    if missing:
        raise ValueError(f"{study}: genes absent from TPM matrix: {', '.join(missing)}")

    records: list[dict[str, object]] = []
    prefix = f"{study}__"
    for column in tpm.columns:
        if not column.startswith(prefix):
            raise ValueError(f"{study}: unexpected TPM column {column}")
        sample_id = column[len(prefix):]
        if sample_id not in metadata.index:
            raise ValueError(f"{study}: metadata missing for {sample_id}")
        meta = metadata.loc[sample_id]
        for gene in GENES:
            value = float(tpm.at[gene, column])
            records.append(
                {
                    "study": study,
                    "sample_id": sample_id,
                    "gene": gene,
                    "condition": str(meta["condition"]),
                    "stage": str(meta.get("stage", "")),
                    "sex": str(meta.get("sex", "")),
                    "tpm": value,
                    "log_tpm": math.log2(value + 1.0),
                }
            )
    return pd.DataFrame.from_records(records)


def load_names(root: Path) -> dict[str, str]:
    table = pd.read_csv(root / "config" / "reference" / "gene_annotations_wbps19.tsv", sep="\t")
    return {
        str(row.gene_id): str(row.functional_name)
        for row in table.itertuples()
        if str(row.gene_id) in GENES
    }


def svg_start(width: int, height: int, title: str, subtitle: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">',
        f"<title>{esc(title)}</title>",
        f"<desc>{esc(subtitle)}</desc>",
        "<style>text{font-family:Inter,Segoe UI,Arial,sans-serif;fill:#172033}.title{font-size:38px;font-weight:750}.panel-title{font-size:23px;font-weight:700}.small{font-size:16px;fill:#536173}.tick{font-size:15px;fill:#536173}.axis{stroke:#94a3b8;stroke-width:1.4}.grid{stroke:#dbe4ee;stroke-width:1.25}.mean{stroke-width:2.8}.rep{fill-opacity:.48;stroke:#fff;stroke-width:1}</style>",
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text class="title" x="80" y="56">{esc(title)}</text>',
    ]


def nice_ceiling(value: float) -> float:
    if value <= 1:
        return 1.0
    return math.ceil(value * 2.0) / 2.0


def y_ticks(ymax: float, count: int = 4) -> list[float]:
    return [i * ymax / count for i in range(count + 1)]


def dot(parts: list[str], x: float, y: float, color: str, shape: str = "circle", radius: float = 4.2, css: str = "rep") -> None:
    if shape == "square":
        parts.append(f'<rect class="{css}" x="{x-radius:.2f}" y="{y-radius:.2f}" width="{2*radius:.2f}" height="{2*radius:.2f}" fill="{color}"/>')
    else:
        parts.append(f'<circle class="{css}" cx="{x:.2f}" cy="{y:.2f}" r="{radius:.2f}" fill="{color}"/>')


def panel_frame(parts: list[str], x: float, y: float, w: float, h: float, target: bool) -> None:
    fill = "#fff7ed" if target else "#fbfdff"
    stroke = "#f97316" if target else "#d9e2ec"
    stroke_width = 2.4 if target else 1.0
    parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}"/>')


def stats(values: pd.Series) -> tuple[float, float]:
    array = values.to_numpy(dtype=float)
    return float(np.mean(array)), float(np.std(array, ddof=1)) if len(array) > 1 else 0.0


def prjeb_figure(data: pd.DataFrame, names: dict[str, str], output: Path) -> None:
    width, height = 2200, 2220
    parts = svg_start(
        width,
        height,
        "PRJEB32839 · expression of six kinases across the complete study",
        "Each category is kept separate. Points are biological samples; large dot and whisker show mean ± SD of log2(TPM + 1).",
    )
    left, right, top, bottom = 190, 80, 88, 275
    gap = 18
    panel_h = (height - top - bottom - gap * 5) / 6
    plot_left = left + 260
    plot_right = width - right - 20
    x_positions = np.linspace(plot_left, plot_right, len(PRJEB_ORDER))
    ymax = nice_ceiling(float(data.log_tpm.max()) * 1.08)
    offsets = np.linspace(-11, 11, 5)

    for row, gene in enumerate(GENES):
        y0 = top + row * (panel_h + gap)
        panel_frame(parts, 70, y0, width - 140, panel_h, gene == TARGET)
        color = "#dc2626" if gene == TARGET else "#2563eb"
        title_color = "#b91c1c" if gene == TARGET else "#172033"
        parts.append(f'<text class="panel-title" x="95" y="{y0+36:.2f}" fill="{title_color}">{esc(gene)}</text>')
        parts.append(f'<text class="small" x="95" y="{y0+62:.2f}">{esc(names.get(gene, "kinase"))}</text>')
        chart_top, chart_bottom = y0 + 22, y0 + panel_h - 38
        for tick in y_ticks(ymax):
            yy = chart_bottom - (tick / ymax) * (chart_bottom - chart_top)
            parts.append(f'<line class="grid" x1="{plot_left}" y1="{yy:.2f}" x2="{plot_right}" y2="{yy:.2f}"/>')
            parts.append(f'<text class="tick" text-anchor="end" x="{plot_left-12}" y="{yy+4:.2f}">{tick:.1f}</text>')
        gene_data = data[data.gene == gene]
        for index, condition in enumerate(PRJEB_ORDER):
            values = gene_data[gene_data.condition == condition].log_tpm
            if values.empty:
                continue
            mean, sd = stats(values)
            x = float(x_positions[index])
            for sample_index, value in enumerate(values):
                jitter = float(offsets[sample_index]) if len(values) == 5 else (sample_index - (len(values)-1)/2) * 5
                yy = chart_bottom - (float(value) / ymax) * (chart_bottom - chart_top)
                dot(parts, x + jitter, yy, color, radius=5.2)
            mean_y = chart_bottom - (mean / ymax) * (chart_bottom - chart_top)
            high_y = chart_bottom - (min(ymax, mean + sd) / ymax) * (chart_bottom - chart_top)
            low_y = chart_bottom - (max(0.0, mean - sd) / ymax) * (chart_bottom - chart_top)
            parts.append(f'<line x1="{x:.2f}" y1="{high_y:.2f}" x2="{x:.2f}" y2="{low_y:.2f}" stroke="{color}" stroke-width="2"/>')
            parts.append(f'<line x1="{x-6:.2f}" y1="{high_y:.2f}" x2="{x+6:.2f}" y2="{high_y:.2f}" stroke="{color}" stroke-width="2"/>')
            parts.append(f'<line x1="{x-6:.2f}" y1="{low_y:.2f}" x2="{x+6:.2f}" y2="{low_y:.2f}" stroke="{color}" stroke-width="2"/>')
            dot(parts, x, mean_y, color, radius=7.7, css="mean")
        parts.append(f'<text class="small" text-anchor="middle" transform="translate(112 {y0+panel_h/2:.2f}) rotate(-90)">log2(TPM + 1)</text>')

    bottom_y = top + 5 * (panel_h + gap) + panel_h
    for x, condition in zip(x_positions, PRJEB_ORDER):
        first, second = LABELS[condition].split("\n")
        parts.append(f'<text class="tick" text-anchor="end" transform="translate({x+5:.2f} {bottom_y+20:.2f}) rotate(-48)"><tspan>{esc(first)}</tspan><tspan x="0" dy="14">{esc(second)}</tspan></text>')
    parts.append(f'<text class="small" x="80" y="{height-28}">Source: frozen HelixForge-SMansoni TPM matrix. Smp_073340 is highlighted in red; no categories were pooled across stage or sex.</text>')
    parts.append("</svg>")
    output.write_text("\n".join(parts), encoding="utf-8")


def trajectory_figure(
    adults: pd.DataFrame,
    names: dict[str, str],
    output: Path,
    juvenile: pd.DataFrame | None = None,
) -> None:
    compare = juvenile is not None
    width, height = 2100, 1450
    title = (
        "Juvenile 26d (PRJEB32839) versus adults (E-ERAD-478) · six kinases"
        if compare
        else "E-ERAD-478 · adult expression trajectories of six kinases"
    )
    subtitle = (
        "Descriptive cross-study TPM view; juvenile sex groups and every adult day/sex/pairing group remain separate."
        if compare
        else "Every day, sex and paired/unpaired group remains separate. Points are samples; trajectories connect exact-condition means."
    )
    parts = svg_start(width, height, title, subtitle)
    panel_w, panel_h = 920, 350
    x_origins = [100, 1080]
    y_origins = [92, 502, 912]
    adult_days = [21, 28, 35, 38]
    x_labels = ["Juvenile\n26d", "Adult\n21d", "Adult\n28d", "Adult\n35d", "Adult\n38d"] if compare else ["21d", "28d", "35d", "38d"]
    all_values = adults.log_tpm
    if compare:
        all_values = pd.concat([all_values, juvenile.log_tpm], ignore_index=True)
    ymax = nice_ceiling(float(all_values.max()) * 1.08)

    for gene_index, gene in enumerate(GENES):
        col, row = gene_index % 2, gene_index // 2
        x0, y0 = x_origins[col], y_origins[row]
        panel_frame(parts, x0, y0, panel_w, panel_h, gene == TARGET)
        title_color = "#b91c1c" if gene == TARGET else "#172033"
        parts.append(f'<text class="panel-title" x="{x0+24}" y="{y0+34}" fill="{title_color}">{esc(gene)} · {esc(names.get(gene, "kinase"))}</text>')
        chart_left, chart_right = x0 + 84, x0 + panel_w - 28
        chart_top, chart_bottom = y0 + 55, y0 + panel_h - 66
        positions = np.linspace(chart_left + 28, chart_right - 28, len(x_labels))
        for tick in y_ticks(ymax):
            yy = chart_bottom - (tick / ymax) * (chart_bottom - chart_top)
            parts.append(f'<line class="grid" x1="{chart_left}" y1="{yy:.2f}" x2="{chart_right}" y2="{yy:.2f}"/>')
            parts.append(f'<text class="tick" text-anchor="end" x="{chart_left-10}" y="{yy+4:.2f}">{tick:.1f}</text>')

        adult_gene = adults[adults.gene == gene]
        adult_offset = 1 if compare else 0
        for series_index, (series, (label, color, shape)) in enumerate(SERIES.items()):
            mean_points: list[tuple[float, float]] = []
            for day_index, day in enumerate(adult_days):
                condition = f"d{day}_{series}"
                values = adult_gene[adult_gene.condition == condition].log_tpm
                if values.empty:
                    continue
                mean, sd = stats(values)
                x = float(positions[day_index + adult_offset])
                spread = np.linspace(-7, 7, len(values))
                for sample_index, value in enumerate(values):
                    yy = chart_bottom - (float(value) / ymax) * (chart_bottom - chart_top)
                    dot(parts, x + float(spread[sample_index]), yy, color, shape=shape, radius=4.7)
                mean_y = chart_bottom - (mean / ymax) * (chart_bottom - chart_top)
                high_y = chart_bottom - (min(ymax, mean + sd) / ymax) * (chart_bottom - chart_top)
                low_y = chart_bottom - (max(0.0, mean - sd) / ymax) * (chart_bottom - chart_top)
                parts.append(f'<line x1="{x:.2f}" y1="{high_y:.2f}" x2="{x:.2f}" y2="{low_y:.2f}" stroke="{color}" stroke-width="1.5"/>')
                mean_points.append((x, mean_y))
            if len(mean_points) > 1:
                points = " ".join(f"{x:.2f},{y:.2f}" for x, y in mean_points)
                parts.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.3"/>')
            for x, y in mean_points:
                dot(parts, x, y, color, shape=shape, radius=7.0, css="mean")

        if compare:
            juvenile_gene = juvenile[juvenile.gene == gene]
            for sex, color, x_shift in [("female", "#e11d48", -9), ("male", "#0f766e", 9)]:
                values = juvenile_gene[juvenile_gene.condition == f"juvenile_26d_{sex}"].log_tpm
                if values.empty:
                    continue
                mean, sd = stats(values)
                x = float(positions[0] + x_shift)
                spread = np.linspace(-5, 5, len(values))
                for sample_index, value in enumerate(values):
                    yy = chart_bottom - (float(value) / ymax) * (chart_bottom - chart_top)
                    dot(parts, x + float(spread[sample_index]), yy, color, shape="circle", radius=4.7)
                mean_y = chart_bottom - (mean / ymax) * (chart_bottom - chart_top)
                high_y = chart_bottom - (min(ymax, mean + sd) / ymax) * (chart_bottom - chart_top)
                low_y = chart_bottom - (max(0.0, mean - sd) / ymax) * (chart_bottom - chart_top)
                parts.append(f'<line x1="{x:.2f}" y1="{high_y:.2f}" x2="{x:.2f}" y2="{low_y:.2f}" stroke="{color}" stroke-width="1.5"/>')
                dot(parts, x, mean_y, color, radius=7.0, css="mean")

        for x, label in zip(positions, x_labels):
            lines = label.split("\n")
            parts.append(f'<text class="tick" text-anchor="middle" x="{x:.2f}" y="{chart_bottom+22:.2f}">{esc(lines[0])}</text>')
            if len(lines) > 1:
                parts.append(f'<text class="tick" text-anchor="middle" x="{x:.2f}" y="{chart_bottom+36:.2f}">{esc(lines[1])}</text>')
        parts.append(f'<text class="small" text-anchor="middle" transform="translate({x0+22} {y0+panel_h/2:.2f}) rotate(-90)">log2(TPM + 1)</text>')

    legend_y = height - 74
    legend_x = 100
    if compare:
        for label, color in [("PRJEB32839 juvenile female", "#e11d48"), ("PRJEB32839 juvenile male", "#0f766e")]:
            dot(parts, legend_x, legend_y, color, radius=6.5, css="mean")
            parts.append(f'<text class="small" x="{legend_x+13}" y="{legend_y+4}">{esc(label)}</text>')
            legend_x += 245
    for _, (label, color, shape) in SERIES.items():
        dot(parts, legend_x, legend_y, color, shape=shape, radius=6.5, css="mean")
        parts.append(f'<text class="small" x="{legend_x+13}" y="{legend_y+4}">{esc(label)}</text>')
        legend_x += 190
    footnote = (
        "Cross-study values are descriptive and were not batch-corrected or tested jointly. Smp_073340 is identified by the orange panel border."
        if compare
        else "Frozen E-ERAD-478 TPM matrix; exact condition means are connected only across adult collection days. Smp_073340 has the orange panel border."
    )
    parts.append(f'<text class="small" x="100" y="{height-28}">{esc(footnote)}</text>')
    parts.append("</svg>")
    output.write_text("\n".join(parts), encoding="utf-8")


def gonad_figure(data: pd.DataFrame, names: dict[str, str], output: Path) -> None:
    """Plot PRJEB14695 gonad samples without pooling tissue or pairing context."""
    data = data[data.condition.isin(GONAD_CONDITIONS)].copy()
    width, height = 2100, 1450
    parts = svg_start(
        width,
        height,
        "PRJEB14695 · expression of six kinases in ovaries and testes",
        "Ovary/testis and mixed-sex/single-sex conditions remain separate. Points are biological samples; large dots show mean ± SD.",
    )
    panel_w, panel_h = 920, 350
    x_origins = [100, 1080]
    y_origins = [92, 502, 912]
    ymax = nice_ceiling(float(data.log_tpm.max()) * 1.08)

    for gene_index, gene in enumerate(GENES):
        col, row = gene_index % 2, gene_index // 2
        x0, y0 = x_origins[col], y_origins[row]
        panel_frame(parts, x0, y0, panel_w, panel_h, gene == TARGET)
        title_color = "#b91c1c" if gene == TARGET else "#172033"
        parts.append(f'<text class="panel-title" x="{x0+24}" y="{y0+34}" fill="{title_color}">{esc(gene)} · {esc(names.get(gene, "kinase"))}</text>')
        chart_left, chart_right = x0 + 84, x0 + panel_w - 28
        chart_top, chart_bottom = y0 + 55, y0 + panel_h - 72
        positions = np.linspace(chart_left + 45, chart_right - 45, len(GONAD_CONDITIONS))
        for tick in y_ticks(ymax):
            yy = chart_bottom - (tick / ymax) * (chart_bottom - chart_top)
            parts.append(f'<line class="grid" x1="{chart_left}" y1="{yy:.2f}" x2="{chart_right}" y2="{yy:.2f}"/>')
            parts.append(f'<text class="tick" text-anchor="end" x="{chart_left-10}" y="{yy+4:.2f}">{tick:.1f}</text>')

        gene_data = data[data.gene == gene]
        for index, condition in enumerate(GONAD_CONDITIONS):
            values = gene_data[gene_data.condition == condition].log_tpm
            if values.empty:
                continue
            mean, sd = stats(values)
            x = float(positions[index])
            color = GONAD_COLORS[condition]
            spread = np.linspace(-9, 9, len(values))
            for sample_index, value in enumerate(values):
                yy = chart_bottom - (float(value) / ymax) * (chart_bottom - chart_top)
                dot(parts, x + float(spread[sample_index]), yy, color, radius=5.2)
            mean_y = chart_bottom - (mean / ymax) * (chart_bottom - chart_top)
            high_y = chart_bottom - (min(ymax, mean + sd) / ymax) * (chart_bottom - chart_top)
            low_y = chart_bottom - (max(0.0, mean - sd) / ymax) * (chart_bottom - chart_top)
            parts.append(f'<line x1="{x:.2f}" y1="{high_y:.2f}" x2="{x:.2f}" y2="{low_y:.2f}" stroke="{color}" stroke-width="2"/>')
            parts.append(f'<line x1="{x-6:.2f}" y1="{high_y:.2f}" x2="{x+6:.2f}" y2="{high_y:.2f}" stroke="{color}" stroke-width="2"/>')
            parts.append(f'<line x1="{x-6:.2f}" y1="{low_y:.2f}" x2="{x+6:.2f}" y2="{low_y:.2f}" stroke="{color}" stroke-width="2"/>')
            dot(parts, x, mean_y, color, radius=7.7, css="mean")

        for x, condition in zip(positions, GONAD_CONDITIONS):
            line1, line2 = GONAD_LABELS[condition].split("\n")
            parts.append(f'<text class="tick" text-anchor="middle" x="{x:.2f}" y="{chart_bottom+24:.2f}">{esc(line1)}</text>')
            parts.append(f'<text class="tick" text-anchor="middle" x="{x:.2f}" y="{chart_bottom+39:.2f}">{esc(line2)}</text>')
        parts.append(f'<text class="small" text-anchor="middle" transform="translate({x0+22} {y0+panel_h/2:.2f}) rotate(-90)">log2(TPM + 1)</text>')

    legend_x, legend_y = 100, height - 74
    for condition in GONAD_CONDITIONS:
        color = GONAD_COLORS[condition]
        label = GONAD_LABELS[condition].replace("\n", " · ")
        dot(parts, legend_x, legend_y, color, radius=6.5, css="mean")
        parts.append(f'<text class="small" x="{legend_x+13}" y="{legend_y+4}">{esc(label)}</text>')
        legend_x += 245
    parts.append(f'<text class="small" x="100" y="{height-28}">Source: frozen PRJEB14695 TPM matrix. Whole-worm samples are excluded; Smp_073340 is identified by the orange panel border.</text>')
    parts.append("</svg>")
    output.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    names = load_names(root)
    prjeb = load_study(root, "PRJEB32839")
    eerad = load_study(root, "E-ERAD-478")
    gonads = load_study(root, "PRJEB14695")

    prjeb_figure(prjeb, names, output / "01_PRJEB32839_complete_six_kinases.svg")
    trajectory_figure(
        eerad,
        names,
        output / "02_juvenile26d_vs_EERAD478_adults_six_kinases.svg",
        juvenile=prjeb[prjeb.condition.isin(["juvenile_26d_female", "juvenile_26d_male"])].copy(),
    )
    trajectory_figure(eerad, names, output / "03_EERAD478_complete_six_kinases.svg")
    gonad_figure(gonads, names, output / "04_PRJEB14695_ovary_testis_six_kinases.svg")

    summary_rows = []
    for frame in (prjeb, eerad, gonads[gonads.condition.isin(GONAD_CONDITIONS)]):
        for (study, gene, condition), group in frame.groupby(["study", "gene", "condition"], sort=False):
            summary_rows.append(
                {
                    "study": study,
                    "gene_id": gene,
                    "condition": condition,
                    "n": len(group),
                    "mean_tpm": group.tpm.mean(),
                    "sd_tpm": group.tpm.std(ddof=1),
                    "mean_log2_tpm_plus_1": group.log_tpm.mean(),
                    "sd_log2_tpm_plus_1": group.log_tpm.std(ddof=1),
                }
            )
    pd.DataFrame(summary_rows).to_csv(output / "kinase_expression_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
