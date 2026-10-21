#!/usr/bin/env python3
"""Cluster differentially expressed epigenetic regulators across RNA-seq studies.

This is an exploratory candidate-prioritization analysis. Differential
expression remains defined by the frozen per-study DESeq2 results. Expression
profiles are standardized within each study before clustering, so raw TPM from
different studies is never treated as a jointly normalized experiment.
"""

from __future__ import annotations

import argparse
import html
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


ALPHA = 0.05
LFC_THRESHOLD = 1.0
PANEL_SIZE = 5
CLUSTER_COLORS = ["#2563eb", "#7c3aed", "#059669", "#d97706", "#db2777"]

ANNOTATION_CONFLICT_TERMS = (
    "mitochondrial ribosomal",
    "ribosomal protein",
    "translation initiation factor",
)
ANNOTATION_SUPPORT_TERMS = (
    "histone",
    "chromobox",
    "heterochromatin",
    "methyltransferase",
    "demethylase",
    "jumonji",
    "jmjc",
    "tudor",
    "chromo domain",
    "phd-type",
    "pwwp",
    "set domain",
    "acetyltransferase",
    "deacetylase",
    "amine oxidase",
)


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def compact(value: float, digits: int = 3) -> str:
    if not math.isfinite(value):
        return "NA"
    return f"{value:.{digits}g}"


def assess_annotation(functional_name: str) -> str:
    value = functional_name.lower()
    if any(term in value for term in ANNOTATION_CONFLICT_TERMS):
        return "CONFLICT"
    if any(term in value for term in ANNOTATION_SUPPORT_TERMS):
        return "SUPPORTED"
    return "REVIEW"


def functional_family(functional_name: str) -> str:
    value = functional_name.lower()
    if "chromobox" in value or "heterochromatin" in value:
        return "Chromobox/HP1"
    if "jumonji" in value or "jmjc" in value or "demethylase" in value or "amine oxidase" in value:
        return "Histone demethylase"
    if "arginine" in value and "methyltransferase" in value:
        return "Arginine methyltransferase"
    if "h4" in value and "20" in value and "methyltransferase" in value:
        return "H4K20 methyltransferase"
    if "setb1" in value:
        return "SETB1-like methyltransferase"
    if "histone-lysine" in value or "set domain" in value or "methyltransferase" in value:
        return "Lysine/SET methyltransferase"
    if "tudor" in value:
        return "Tudor reader"
    if "pwwp" in value:
        return "PWWP reader"
    if "phd-type" in value:
        return "PHD reader"
    if "acetyltransferase" in value:
        return "Histone acetyltransferase"
    if "deacetylase" in value:
        return "Histone deacetylase"
    return "Other/uncertain"


def interpolate(start: tuple[int, int, int], end: tuple[int, int, int], fraction: float) -> str:
    fraction = max(0.0, min(1.0, fraction))
    rgb = tuple(round(a + (b - a) * fraction) for a, b in zip(start, end))
    return "#" + "".join(f"{channel:02x}" for channel in rgb)


def diverging_color(value: float, limit: float) -> str:
    if not math.isfinite(value):
        return "#e5e7eb"
    value = max(-limit, min(limit, value))
    white = (255, 255, 255)
    blue = (49, 54, 149)
    red = (165, 0, 38)
    return interpolate(blue, white, (value + limit) / limit) if value <= 0 else interpolate(white, red, value / limit)


@dataclass
class Node:
    members: tuple[int, ...]
    height: float
    left: int | None = None
    right: int | None = None


def average_linkage(distance: np.ndarray) -> tuple[dict[int, Node], int]:
    n = distance.shape[0]
    nodes: dict[int, Node] = {index: Node((index,), 0.0) for index in range(n)}
    active = list(range(n))
    next_id = n
    while len(active) > 1:
        best: tuple[float, int, int] | None = None
        for offset, left in enumerate(active):
            for right in active[offset + 1 :]:
                values = distance[np.ix_(nodes[left].members, nodes[right].members)]
                candidate = (float(np.mean(values)), left, right)
                if best is None or candidate < best:
                    best = candidate
        assert best is not None
        height, left, right = best
        nodes[next_id] = Node(nodes[left].members + nodes[right].members, height, left, right)
        active = [item for item in active if item not in {left, right}] + [next_id]
        next_id += 1
    return nodes, active[0]


def ordered_leaves(nodes: dict[int, Node], node_id: int, distance: np.ndarray) -> list[int]:
    node = nodes[node_id]
    if node.left is None or node.right is None:
        return [node.members[0]]
    left = ordered_leaves(nodes, node.left, distance)
    right = ordered_leaves(nodes, node.right, distance)
    candidates = []
    for first, second in ((left, right), (right, left)):
        for first_order in (first, list(reversed(first))):
            for second_order in (second, list(reversed(second))):
                boundary = distance[first_order[-1], second_order[0]]
                candidates.append((float(boundary), first_order + second_order))
    return min(candidates, key=lambda item: (item[0], item[1]))[1]


def cut_tree(nodes: dict[int, Node], root: int, k: int, leaf_order: list[int]) -> tuple[dict[int, int], list[int]]:
    active = [root]
    while len(active) < k:
        splittable = [item for item in active if nodes[item].left is not None]
        if not splittable:
            break
        split = max(splittable, key=lambda item: nodes[item].height)
        active.remove(split)
        active.extend([nodes[split].left, nodes[split].right])
    position = {leaf: index for index, leaf in enumerate(leaf_order)}
    active.sort(key=lambda item: min(position[member] for member in nodes[item].members))
    labels: dict[int, int] = {}
    for cluster_id, item in enumerate(active, start=1):
        for member in nodes[item].members:
            labels[member] = cluster_id
    return labels, active


def silhouette(distance: np.ndarray, labels: dict[int, int]) -> float:
    values = []
    n = distance.shape[0]
    for index in range(n):
        own = labels[index]
        own_members = [item for item in range(n) if labels[item] == own and item != index]
        if not own_members:
            values.append(0.0)
            continue
        a = float(np.mean(distance[index, own_members]))
        alternatives = []
        for cluster in sorted(set(labels.values()) - {own}):
            members = [item for item in range(n) if labels[item] == cluster]
            alternatives.append(float(np.mean(distance[index, members])))
        b = min(alternatives)
        values.append((b - a) / max(a, b) if max(a, b) else 0.0)
    return float(np.mean(values))


def parse_candidate_metadata(selected: pd.DataFrame) -> pd.DataFrame:
    records: dict[str, dict[str, object]] = {}
    for row in selected.itertuples(index=False):
        resolved = str(row.resolved_gene_id)
        if not resolved:
            continue
        record = records.setdefault(
            resolved,
            {
                "gene_id": resolved,
                "requested_gene_ids": [],
                "groups": set(),
                "functional_name": str(row.functional_name),
                "display_label": str(row.display_label),
                "resolution_statuses": set(),
            },
        )
        record["requested_gene_ids"].append(str(row.gene_id))
        record["groups"].update(item for item in str(row.groups).split(";") if item)
        record["resolution_statuses"].add(str(row.resolution_status))
    output = []
    for record in records.values():
        output.append(
            {
                **record,
                "requested_gene_ids": ";".join(dict.fromkeys(record["requested_gene_ids"])),
                "groups": ";".join(sorted(record["groups"])),
                "resolution_statuses": ";".join(sorted(record["resolution_statuses"])),
            }
        )
    return pd.DataFrame(output).set_index("gene_id")


def load_expression_profiles(
    root: Path,
    studies: list[str],
    genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw_blocks = []
    z_blocks = []
    variability_records = []
    for study in studies:
        tpm = pd.read_csv(root / "results" / study / "expression" / "tpm_matrix.tsv", sep="\t", index_col=0)
        metadata = pd.read_csv(root / "metadata" / study / "samples.tsv", sep="\t")
        available = [gene for gene in genes if gene in tpm.index]
        if len(available) != len(genes):
            missing = sorted(set(genes) - set(available))
            raise ValueError(f"{study}: missing clustered genes: {missing}")
        metadata = metadata.set_index("sample_id")
        conditions = list(dict.fromkeys(str(value) for value in metadata.condition))
        condition_means = []
        for condition in conditions:
            sample_ids = [str(sample) for sample in metadata[metadata.condition == condition].index]
            columns = [f"{study}__{sample}" for sample in sample_ids]
            absent = [column for column in columns if column not in tpm.columns]
            if absent:
                raise ValueError(f"{study}: TPM columns absent: {absent}")
            values = tpm.loc[genes, columns].astype(float)
            log_values = np.log2(values + 1.0)
            condition_means.append(values.mean(axis=1).rename(f"{study}:{condition}"))
            for gene in genes:
                gene_values = values.loc[gene].to_numpy(dtype=float)
                mean = float(np.mean(gene_values))
                variability_records.append(
                    {
                        "gene_id": gene,
                        "context": f"{study}:{condition}",
                        "study": study,
                        "condition": condition,
                        "n": len(gene_values),
                        "mean_tpm": mean,
                        "sd_tpm": float(np.std(gene_values, ddof=1)) if len(gene_values) > 1 else 0.0,
                        "cv_tpm": float(np.std(gene_values, ddof=1) / mean) if len(gene_values) > 1 and mean > 0 else math.nan,
                        "mean_log2_tpm_plus_1": float(log_values.loc[gene].mean()),
                    }
                )
        raw_block = pd.concat(condition_means, axis=1)
        log_block = np.log2(raw_block + 1.0)
        centered = log_block.sub(log_block.mean(axis=1), axis=0)
        sd = log_block.std(axis=1, ddof=0).replace(0.0, 1.0)
        z_block = centered.div(sd, axis=0) / math.sqrt(len(conditions))
        raw_blocks.append(raw_block)
        z_blocks.append(z_block)
    return pd.concat(raw_blocks, axis=1), pd.concat(z_blocks, axis=1), pd.DataFrame(variability_records)


def dendrogram_segments(nodes: dict[int, Node], root: int, order: list[int], left: float, right: float, y_positions: dict[int, float]) -> list[str]:
    parts: list[str] = []
    max_height = max(nodes[root].height, 1e-9)

    def visit(node_id: int) -> tuple[float, float]:
        node = nodes[node_id]
        if node.left is None or node.right is None:
            return right, y_positions[node.members[0]]
        left_x, left_y = visit(node.left)
        right_x, right_y = visit(node.right)
        x = right - (node.height / max_height) * (right - left)
        parts.append(f'<line x1="{x:.2f}" y1="{left_y:.2f}" x2="{left_x:.2f}" y2="{left_y:.2f}" stroke="#64748b" stroke-width="1.2"/>')
        parts.append(f'<line x1="{x:.2f}" y1="{right_y:.2f}" x2="{right_x:.2f}" y2="{right_y:.2f}" stroke="#64748b" stroke-width="1.2"/>')
        parts.append(f'<line x1="{x:.2f}" y1="{left_y:.2f}" x2="{x:.2f}" y2="{right_y:.2f}" stroke="#64748b" stroke-width="1.2"/>')
        return x, (left_y + right_y) / 2.0

    visit(root)
    return parts


def heatmap_svg(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    row_order: list[str],
    cluster_by_gene: dict[str, int],
    nodes: dict[int, Node],
    root_node: int,
    index_by_gene: dict[str, int],
    title: str,
    legend_label: str,
    value_limit: float,
    cluster_count: int,
    significance: pd.DataFrame | None = None,
) -> str:
    rows = len(row_order)
    columns = len(matrix.columns)
    cell_w = 31 if columns > 40 else 55
    cell_h = 31
    dendro_left = 70
    dendro_right = 285
    label_left = 305
    heat_left = 790
    heat_top = 115
    width = heat_left + columns * cell_w + 130
    height = heat_top + rows * cell_h + 480
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">',
        f'<title>{esc(title)}</title>',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:Inter,Segoe UI,Arial,sans-serif;fill:#172033}.title{font-size:34px;font-weight:750}.row{font-size:14px}.small{font-size:13px;fill:#536173}.col{font-size:12px;fill:#334155}</style>',
        f'<text class="title" x="70" y="50">{esc(title)}</text>',
    ]
    y_positions = {index_by_gene[gene]: heat_top + row_index * cell_h + cell_h / 2 for row_index, gene in enumerate(row_order)}
    parts.extend(dendrogram_segments(nodes, root_node, [index_by_gene[g] for g in row_order], dendro_left, dendro_right, y_positions))

    for row_index, gene in enumerate(row_order):
        y = heat_top + row_index * cell_h
        cluster = cluster_by_gene[gene]
        cluster_color = CLUSTER_COLORS[(cluster - 1) % len(CLUSTER_COLORS)]
        parts.append(f'<rect x="{label_left}" y="{y+2}" width="10" height="{cell_h-4}" fill="{cluster_color}"/>')
        functional = str(metadata.at[gene, "functional_name"])
        groups = str(metadata.at[gene, "groups"])
        label = f"{gene} · {functional} [{groups}]"
        if len(label) > 66:
            label = label[:63].rstrip() + "…"
        parts.append(f'<text class="row" x="{label_left+18}" y="{y+20}">{esc(label)}</text>')
        for column_index, column in enumerate(matrix.columns):
            value = float(matrix.at[gene, column])
            x = heat_left + column_index * cell_w
            parts.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{diverging_color(value, value_limit)}" stroke="#ffffff" stroke-width="0.7"><title>{esc(gene)} | {esc(column)} | {compact(value, 4)}</title></rect>')
            if significance is not None and bool(significance.at[gene, column]):
                parts.append(f'<circle cx="{x+cell_w/2:.2f}" cy="{y+cell_h/2:.2f}" r="4.2" fill="#111827"/>')
    bottom = heat_top + rows * cell_h
    for column_index, column in enumerate(matrix.columns):
        x = heat_left + column_index * cell_w + cell_w / 2
        parts.append(f'<text class="col" text-anchor="end" transform="translate({x:.2f} {bottom+15}) rotate(-58)">{esc(column)}</text>')

    legend_y = height - 70
    for index in range(101):
        value = -value_limit + (2 * value_limit * index / 100)
        parts.append(f'<rect x="{heat_left + index*3}" y="{legend_y}" width="3" height="18" fill="{diverging_color(value, value_limit)}"/>')
    parts.append(f'<text class="small" x="{heat_left}" y="{legend_y+38}">−{value_limit:g}</text>')
    parts.append(f'<text class="small" text-anchor="end" x="{heat_left+303}" y="{legend_y+38}">{value_limit:g}</text>')
    parts.append(f'<text class="small" x="{heat_left+330}" y="{legend_y+16}">{esc(legend_label)}</text>')
    if significance is not None:
        parts.append(f'<circle cx="{heat_left+650}" cy="{legend_y+9}" r="4.2" fill="#111827"/>')
        parts.append(f'<text class="small" x="{heat_left+663}" y="{legend_y+14}">padj &lt; 0.05 and |log2FC| ≥ 1</text>')
    for cluster in range(1, cluster_count + 1):
        x = 70 + (cluster - 1) * 125
        color = CLUSTER_COLORS[cluster - 1]
        parts.append(f'<rect x="{x}" y="{legend_y}" width="14" height="14" fill="{color}"/>')
        parts.append(f'<text class="small" x="{x+20}" y="{legend_y+13}">Cluster {cluster}</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def ranking_svg(ranking: pd.DataFrame, recommended: set[str]) -> str:
    shown = ranking.sort_values("priority_score", ascending=False).head(15).iloc[::-1]
    width, height = 1900, 1040
    left, right, top, bottom = 650, 110, 95, 100
    plot_w = width - left - right
    row_h = (height - top - bottom) / len(shown)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">',
        '<title>Epigenetic regulator prioritization</title>',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:Inter,Segoe UI,Arial,sans-serif;fill:#172033}.title{font-size:36px;font-weight:750}.label{font-size:17px}.score{font-size:15px;font-weight:700}.small{font-size:14px;fill:#536173}</style>',
        '<text class="title" x="70" y="52">Exploratory prioritization of differentially expressed epigenetic regulators</text>',
    ]
    for tick in range(0, 101, 20):
        x = left + plot_w * tick / 100
        parts.append(f'<line x1="{x}" y1="{top}" x2="{x}" y2="{height-bottom}" stroke="#e2e8f0"/>')
        parts.append(f'<text class="small" text-anchor="middle" x="{x}" y="{height-bottom+28}">{tick}</text>')
    for index, row in enumerate(shown.itertuples()):
        y = top + index * row_h
        color = CLUSTER_COLORS[(int(row.cluster) - 1) % len(CLUSTER_COLORS)]
        bar_w = plot_w * float(row.priority_score) / 100
        opacity = 1.0 if row.gene_id in recommended else 0.62
        parts.append(f'<rect x="{left}" y="{y+8:.2f}" width="{bar_w:.2f}" height="{row_h-16:.2f}" rx="5" fill="{color}" opacity="{opacity}"/>')
        star = "★ " if row.gene_id in recommended else ""
        label = f"{star}{row.gene_id} · {row.functional_name} [{row.groups}]"
        if len(label) > 72:
            label = label[:69].rstrip() + "…"
        weight = "700" if row.gene_id in recommended else "500"
        parts.append(f'<text class="label" x="{left-18}" y="{y+row_h/2+6:.2f}" text-anchor="end" font-weight="{weight}">{esc(label)}</text>')
        parts.append(f'<text class="score" x="{left+bar_w+10:.2f}" y="{y+row_h/2+5:.2f}">{row.priority_score:.1f}</text>')
    parts.append(f'<text class="small" x="70" y="{height-30}">★ preliminary five-gene panel: strict DE, coherent annotation and functional diversity. Score is a transparent prioritization heuristic, not a p-value.</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    atlas_data = root / "results" / "atlas" / "rnaseq" / "data"
    config = json.loads((root / "analyses" / "atlas" / "atlas_config.json").read_text(encoding="utf-8"))
    studies = [str(item["id"]) for item in config["studies"]]
    selected = pd.read_csv(atlas_data / "selected_genes.tsv", sep="\t", keep_default_na=False)
    candidate_meta = parse_candidate_metadata(selected)
    lfc = pd.read_csv(atlas_data / "gene_contrast_log2fc.tsv", sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce")
    padj = pd.read_csv(atlas_data / "gene_contrast_padj.tsv", sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce")
    catalog = pd.read_csv(atlas_data / "contrast_catalog.tsv", sep="\t", keep_default_na=False)
    catalog = catalog[catalog.contrast_key.isin(lfc.columns)].copy()
    catalog = catalog.set_index("contrast_key", drop=False)

    resolved = [gene for gene in candidate_meta.index if gene in lfc.index]
    strict = (padj.loc[resolved] < ALPHA) & (lfc.loc[resolved].abs() >= LFC_THRESHOLD)
    de_genes = [gene for gene in resolved if bool(strict.loc[gene].any())]
    if len(de_genes) < PANEL_SIZE:
        raise ValueError("too few strict DE epigenetic candidates for five-gene selection")

    raw_expression, expression_features, variability = load_expression_profiles(root, studies, de_genes)
    features = expression_features.loc[de_genes].to_numpy(dtype=float)
    correlation = np.corrcoef(features)
    correlation = np.nan_to_num(correlation, nan=0.0)
    distance = np.clip(1.0 - correlation, 0.0, 2.0)
    np.fill_diagonal(distance, 0.0)
    nodes, root_node = average_linkage(distance)
    leaf_indices = ordered_leaves(nodes, root_node, distance)
    leaf_genes = [de_genes[index] for index in leaf_indices]
    index_by_gene = {gene: index for index, gene in enumerate(de_genes)}

    silhouette_rows = []
    for k in range(2, min(9, len(de_genes))):
        trial_labels, _ = cut_tree(nodes, root_node, k, leaf_indices)
        silhouette_rows.append({"clusters": k, "mean_silhouette": silhouette(distance, trial_labels)})
    best_k = int(max(silhouette_rows, key=lambda row: row["mean_silhouette"])["clusters"])
    labels, _ = cut_tree(nodes, root_node, best_k, leaf_indices)
    cluster_by_gene = {de_genes[index]: labels[index] for index in range(len(de_genes))}
    pd.DataFrame(silhouette_rows).to_csv(output / "cluster_silhouette_diagnostic.tsv", sep="\t", index=False)

    contrast_study = catalog.study.to_dict()
    ranking_rows = []
    for gene in de_genes:
        gene_sig = strict.loc[gene]
        sig_columns = [column for column, value in gene_sig.items() if bool(value)]
        sig_studies = sorted(set(str(contrast_study[column]) for column in sig_columns))
        sig_effects = lfc.loc[gene, sig_columns].abs().dropna()
        context_values = raw_expression.loc[gene]
        gene_variability = variability[variability.gene_id == gene]
        median_cv = float(gene_variability.cv_tpm.median(skipna=True))
        fraction_detected = float(np.mean(context_values.to_numpy(dtype=float) > 1.0))
        median_tpm = float(np.median(context_values.to_numpy(dtype=float)))
        max_abs_lfc = float(sig_effects.max())
        functional = str(candidate_meta.at[gene, "functional_name"])
        annotation_status = assess_annotation(functional)
        family = functional_family(functional)
        annotation_specific = bool(functional and functional != gene and "uncharacterized" not in functional.lower() and "hypothetical" not in functional.lower())
        score = (
            35.0 * min(len(sig_studies) / 3.0, 1.0)
            + 20.0 * min(math.log1p(len(sig_columns)) / math.log(9.0), 1.0)
            + 20.0 * min(max_abs_lfc / 4.0, 1.0)
            + 10.0 * min(median_tpm / 20.0, 1.0)
            + 7.0 * fraction_detected
            + 5.0 * max(0.0, 1.0 - min(median_cv if math.isfinite(median_cv) else 1.0, 1.0))
            + 3.0 * float(annotation_specific)
        )
        top_contexts = context_values.sort_values(ascending=False).head(3)
        strongest = lfc.loc[gene, sig_columns].abs().sort_values(ascending=False).head(3).index
        ranking_rows.append(
            {
                "gene_id": gene,
                "requested_gene_ids": candidate_meta.at[gene, "requested_gene_ids"],
                "groups": candidate_meta.at[gene, "groups"],
                "functional_name": functional,
                "annotation_status": annotation_status,
                "functional_family": family,
                "cluster": cluster_by_gene[gene],
                "priority_score": score,
                "strict_significant_studies": len(sig_studies),
                "significant_study_ids": ";".join(sig_studies),
                "strict_significant_contrasts": len(sig_columns),
                "max_abs_log2fc": max_abs_lfc,
                "median_context_tpm": median_tpm,
                "fraction_contexts_tpm_gt_1": fraction_detected,
                "median_within_context_cv": median_cv,
                "top_expression_contexts": ";".join(f"{key}={compact(value,4)}" for key, value in top_contexts.items()),
                "strongest_significant_contrasts": ";".join(f"{key}={compact(float(lfc.at[gene, key]),4)}" for key in strongest),
            }
        )
    ranking = pd.DataFrame(ranking_rows).sort_values(["priority_score", "gene_id"], ascending=[False, True])
    eligible = ranking[ranking.annotation_status == "SUPPORTED"].copy()
    selected_ids: list[str] = []
    selection_reasons: dict[str, str] = {}

    def choose(candidates: pd.DataFrame, reason: str, prefer_novel_family: bool = True) -> None:
        if len(selected_ids) >= PANEL_SIZE or candidates.empty:
            return
        candidates = candidates[~candidates.gene_id.isin(selected_ids)].copy()
        if candidates.empty:
            return
        used_families = set(ranking[ranking.gene_id.isin(selected_ids)].functional_family)
        candidates["novel_family"] = ~candidates.functional_family.isin(used_families)
        sort_columns = ["priority_score", "gene_id"]
        ascending = [False, True]
        if prefer_novel_family:
            sort_columns = ["novel_family"] + sort_columns
            ascending = [False] + ascending
        chosen = candidates.sort_values(sort_columns, ascending=ascending).iloc[0]
        selected_ids.append(str(chosen.gene_id))
        selection_reasons[str(chosen.gene_id)] = reason

    # Cover the mechanistic groups first. A multi-group gene can satisfy more than
    # one requirement, avoiding artificial duplication of closely related targets.
    for group in ("Writers", "Erasers", "Grupo HP1", "Readers"):
        already_covered = any(group in str(ranking[ranking.gene_id == gene].iloc[0].groups).split(";") for gene in selected_ids)
        if already_covered:
            continue
        candidates = eligible[eligible.groups.map(lambda value: group in str(value).split(";"))]
        choose(candidates, f"melhor candidato anotado cobrindo {group}")

    # Preserve each natural expression cluster supported by the silhouette.
    for cluster in range(1, best_k + 1):
        if any(int(ranking[ranking.gene_id == gene].iloc[0].cluster) == cluster for gene in selected_ids):
            continue
        choose(eligible[eligible.cluster == cluster], f"melhor representante anotado do cluster de expressão {cluster}")

    # Complete the panel with the strongest candidates from still-unrepresented
    # molecular families whenever possible.
    while len(selected_ids) < PANEL_SIZE:
        before = len(selected_ids)
        choose(eligible, "melhor candidato restante de uma família funcional distinta")
        if len(selected_ids) == before:
            break
    if len(selected_ids) < PANEL_SIZE:
        raise ValueError("fewer than five annotation-supported strict-DE candidates")

    ranking["selection_reason"] = ranking.gene_id.map(selection_reasons).fillna("")
    recommended = ranking[ranking.gene_id.isin(selected_ids)].copy()
    recommended["selection_order"] = recommended.gene_id.map({gene: index + 1 for index, gene in enumerate(selected_ids)})
    recommended = recommended.sort_values("selection_order")
    recommended_ids = set(recommended.gene_id)
    ranking["recommended_five"] = ranking.gene_id.isin(recommended_ids).map({True: "YES", False: "NO"})

    ranking.to_csv(output / "epigenetic_de_candidate_ranking.tsv", sep="\t", index=False)
    recommended.to_csv(output / "recommended_five_genes.tsv", sep="\t", index=False)
    raw_expression.loc[leaf_genes].to_csv(output / "condition_mean_tpm.tsv", sep="\t", index_label="gene_id")
    expression_features.loc[leaf_genes].to_csv(output / "within_study_expression_zscore_weighted.tsv", sep="\t", index_label="gene_id")
    lfc.loc[leaf_genes].to_csv(output / "all_contrast_log2fc.tsv", sep="\t", index_label="gene_id", na_rep="NA")
    padj.loc[leaf_genes].to_csv(output / "all_contrast_padj.tsv", sep="\t", index_label="gene_id", na_rep="NA")
    variability.to_csv(output / "within_condition_expression_variability.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"gene_id": gene, "cluster": cluster_by_gene[gene], "display_order": index + 1} for index, gene in enumerate(leaf_genes)]
    ).to_csv(output / "cluster_membership.tsv", sep="\t", index=False)

    selected_contrasts = []
    for study in list(dict.fromkeys(catalog.study)):
        study_columns = [column for column in catalog[catalog.study == study].contrast_key if column in lfc.columns]
        scored = []
        for position, column in enumerate(study_columns):
            n_sig = int(strict.loc[de_genes, column].sum())
            variance = float(lfc.loc[de_genes, column].var(skipna=True))
            scored.append((n_sig, variance, -position, column))
        chosen = {item[3] for item in sorted(scored, reverse=True)[:4]}
        selected_contrasts.extend(column for column in study_columns if column in chosen)

    expression_svg = heatmap_svg(
        expression_features.loc[leaf_genes],
        candidate_meta,
        leaf_genes,
        cluster_by_gene,
        nodes,
        root_node,
        index_by_gene,
        "Expression clusters of strict-DE epigenetic regulators",
        "within-study z-score (study-balanced)",
        2.5,
        best_k,
    )
    effect_svg = heatmap_svg(
        lfc.loc[leaf_genes, selected_contrasts].fillna(0.0),
        candidate_meta,
        leaf_genes,
        cluster_by_gene,
        nodes,
        root_node,
        index_by_gene,
        "Differential-effect profiles of clustered epigenetic regulators",
        "DESeq2 log2 fold change",
        4.0,
        best_k,
        strict.loc[leaf_genes, selected_contrasts],
    )
    (output / "01_expression_clusters.svg").write_text(expression_svg, encoding="utf-8")
    (output / "02_differential_effect_clusters.svg").write_text(effect_svg, encoding="utf-8")
    (output / "03_candidate_prioritization.svg").write_text(ranking_svg(ranking, recommended_ids), encoding="utf-8")

    unresolved = selected[selected.resolution_status.isin(["AMBIGUOUS_PREVIOUS_STABLE_ID", "NOT_FOUND_WBPS19"])]
    best_silhouette = max(silhouette_rows, key=lambda row: row["mean_silhouette"])
    lines = [
        "# Clustering exploratório de reguladores epigenéticos diferencialmente expressos",
        "",
        "## Resumo",
        "",
        f"- Genes solicitados: **{len(selected)}**.",
        f"- IDs resolvidos de forma não ambígua no WBPS19: **{len(candidate_meta)}**.",
        f"- Reguladores com pelo menos um contraste atendendo `padj < {ALPHA}` e `|log2FC| ≥ {LFC_THRESHOLD}`: **{len(de_genes)}**.",
        f"- Estrutura de clusters sustentada pelo melhor silhouette exploratório: **k={best_k}** (`{best_silhouette['mean_silhouette']:.3f}`).",
        "- O painel de cinco genes não corresponde a cinco clusters forçados. A seleção combina evidência de DE, grupos funcionais, famílias moleculares distintas e coerência da anotação.",
        "",
        "O clustering usa médias de `log2(TPM + 1)` padronizadas separadamente dentro de cada estudo. Cada bloco de estudo recebe peso proporcional a `1/sqrt(n de condições)`, impedindo que estudos com mais condições dominem a distância. A significância continua vindo exclusivamente dos modelos DESeq2 por estudo.",
        "",
        "## Painel preliminar de cinco genes",
        "",
        "| Gene atual | ID da lista original | Função anotada | Grupo | Família | Cluster | Estudos DE | Contrastes DE | Máx. abs(log2FC) | TPM mediano | Escore | Motivo |",
        "|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in recommended.itertuples():
        lines.append(
            f"| `{row.gene_id}` | `{row.requested_gene_ids}` | {row.functional_name} | {row.groups} | {row.functional_family} | {int(row.cluster)} | {int(row.strict_significant_studies)} | {int(row.strict_significant_contrasts)} | {row.max_abs_log2fc:.2f} | {row.median_context_tpm:.2f} | {row.priority_score:.1f} | {row.selection_reason} |"
        )
    lines.extend(
        [
            "",
            "Estes cinco genes são uma **pré-seleção analítica**. O escore combina suporte em estudos e contrastes, magnitude do efeito, expressão detectável, variabilidade entre réplicas e qualidade mínima da anotação. A composição final também exige anotação compatível e diversidade funcional. O escore não é um p-valor.",
            "",
            "Genes com anotação `REVIEW` ou `CONFLICT` permanecem no ranking e nas figuras, mas não entram automaticamente no painel. Isso evita promover um outlier ou um ID legado biologicamente incompatível apenas por possuir escore alto.",
            "",
            "Antes da validação de bancada ainda devem ser verificados: especificidade dos primers, isoformas, unicidade no genoma, tamanho do amplicon, eficiência de amplificação e adequação ao tecido/estágio que será coletado.",
            "",
            "## Exclusões de identificação",
            "",
        ]
    )
    for row in unresolved.itertuples(index=False):
        lines.append(f"- `{row.gene_id}`: `{row.resolution_status}`; excluído do clustering para não atribuir evidência ao gene errado.")
    lines.extend(
        [
            "",
            "## Arquivos principais",
            "",
            "- `01_expression_clusters.svg/png`: agrupamento de padrões de expressão.",
            "- `02_differential_effect_clusters.svg/png`: efeitos DESeq2; pontos pretos marcam o critério estrito.",
            "- `03_candidate_prioritization.svg/png`: ranking transparente; estrelas marcam o painel preliminar.",
            "- `recommended_five_genes.tsv`: painel e métricas completas.",
            "- `epigenetic_de_candidate_ranking.tsv`: todos os candidatos elegíveis.",
            "- Matrizes TSV: dados numéricos reproduzíveis usados nas figuras.",
            "",
            "## Limitações",
            "",
            "- Não foi ajustado um modelo conjunto entre estudos.",
            "- Contrastes representam perguntas biológicas diferentes; o clustering descreve padrões, não causalidade.",
            "- A escolha definitiva depende do contexto experimental da validação.",
        ]
    )
    (output / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    manifest = {
        "status": "PASS",
        "analysis": "epigenetic_strict_de_expression_clustering",
        "alpha": ALPHA,
        "absolute_log2_fold_change": LFC_THRESHOLD,
        "requested_genes": len(selected),
        "resolved_unique_genes": len(candidate_meta),
        "strict_de_genes": len(de_genes),
        "natural_expression_clusters": best_k,
        "selection_strategy": "strict_de_plus_annotation_and_functional_diversity",
        "recommended_genes": list(recommended.gene_id),
        "cross_study_joint_model": False,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
