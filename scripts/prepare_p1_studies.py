#!/usr/bin/env python3
"""Build the pending P1 RNA-seq input packages from frozen ENA snapshots."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PREPARED_ON = "2026-10-09"
ENA_FIELDS = (
    "study_accession,secondary_study_accession,run_accession,experiment_accession,"
    "sample_accession,secondary_sample_accession,sample_alias,sample_title,"
    "sample_description,experiment_title,library_name,library_layout,library_strategy,"
    "library_source,library_selection,instrument_platform,instrument_model,fastq_ftp,"
    "fastq_md5,fastq_bytes,read_count,base_count,first_public,tax_id,scientific_name"
)

STUDIES = {
    "PRJNA312093": {
        "title": "Picard et al. — diferenciação sexual ao longo do desenvolvimento",
        "paper": "https://doi.org/10.1371/journal.pntd.0004930",
        "mode": "full",
        "layout": "SINGLE",
        "expected_runs": 20,
        "expected_samples": 20,
        "fragment_mean": "200",
        "fragment_sd": "80",
        "fragment_source": "prior provisório baseado no protocolo; validar por sensibilidade no piloto biológico",
        "limitations": "O depósito reúne 20 bibliotecas RNA-seq e 24 ChIP-seq sob um único BioSample. O desenho foi reconstruído pelos nomes de biblioteca. O prior de fragmento não foi publicado e permanece um gate biológico.",
    },
    "PRJNA494090": {
        "title": "sensibilidade reduzida ao praziquantel por seleção laboratorial",
        "paper": "https://doi.org/10.1016/j.molbiopara.2018.12.005",
        "mode": "full",
        "layout": "PAIRED",
        "expected_runs": 24,
        "expected_samples": 24,
        "limitations": "Os contrastes comparam a população selecionada e o controle dentro de cada passagem; não representam exposição aguda ao PZQ.",
    },
    "PRJNA704646": {
        "title": "base genética da resposta ao praziquantel",
        "paper": "https://doi.org/10.1126/scitranslmed.abj9114",
        "mode": "import",
        "layout": "PAIRED",
        "expected_runs": 24,
        "expected_samples": 8,
        "limitations": "As três bibliotecas de cada combinação foram produzidas do mesmo pool de RNA. Elas são réplicas técnicas e serão consolidadas; depois da consolidação resta n=1 biológico por combinação, portanto DE inferencial fica bloqueada.",
    },
    "PRJNA1150591": {
        "title": "meclonazepam em juvenis e adultos",
        "paper": "https://doi.org/10.1371/journal.pntd.0012969",
        "mode": "full",
        "layout": "PAIRED",
        "expected_runs": 15,
        "expected_samples": 15,
        "limitations": "O depósito não contém controle de 4 semanas; o artigo informa que esse grupo deve ser pareado à coorte de 4 semanas já processada em PRJNA597909. Nesta execução isolada, o único contraste inferencial pré-registrado é MCLZ versus controle em adultos de 7 semanas; a comparação de 4 semanas fica para a integração posterior entre depósitos.",
    },
    "PRJNA1177852": {
        "title": "curso temporal de meclonazepam",
        "paper": "https://doi.org/10.1371/journal.pntd.0012969",
        "mode": "import",
        "layout": "PAIRED",
        "expected_runs": 8,
        "expected_samples": 8,
        "limitations": "Há uma biblioteca por tempo e nenhum conjunto paralelo de controles. O depósito será quantificado como série descritiva e ligado ao estudo publicado, sem DE inferencial isolada.",
    },
    "PRJNA551487": {
        "title": "resposta ao agonista de canal de cálcio FPL-64176",
        "paper": "https://doi.org/10.1016/j.ijpddr.2019.08.006",
        "mode": "full",
        "layout": "PAIRED",
        "expected_runs": 10,
        "expected_samples": 10,
        "limitations": "A inclusão quantitativa continua condicionada ao gate de QC independente; se falhar, o estudo permanece apenas como evidência publicada.",
    },
    "PRJNA428470": {
        "title": "perturbação por 5-azacitidina em fêmeas adultas",
        "paper": "https://doi.org/10.1038/s41598-020-78669-5",
        "mode": "full",
        "layout": "PAIRED",
        "expected_runs": 6,
        "expected_samples": 6,
        "limitations": "O contraste 3 × 3 testa uma perturbação epigenética; não deve ser interpretado como prova direta de mecanismo de modificação de histonas.",
    },
}

RUN_FIELDS = [
    "dataset", "sample_id", "run_accession", "biosample", "secondary_biosample",
    "experiment_accession", "technical_unit", "library_name", "library_layout",
    "fastq_1_url", "fastq_2_url", "fastq_1_md5", "fastq_2_md5",
    "fastq_1_bytes", "fastq_2_bytes", "read_count", "base_count",
    "ena_first_public", "ena_query_url",
]
SAMPLE_FIELDS = [
    "dataset", "sample_id", "source_sample_name", "biosample", "condition", "stage",
    "tissue", "sex", "infection_mode", "treatment", "time_hours", "batch",
    "replicate", "technical_runs",
]
METADATA_FIELDS = [
    "dataset", "sample_id", "file_prefix", "run_accession", "condition", "stage",
    "tissue", "sex", "infection_mode", "treatment", "time_hours", "batch",
    "replicate", "biosample", "library_name", "technical_unit", "library_layout",
    "fragment_length_mean", "fragment_length_sd",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_rows(path: Path, rows: list[dict[str, str]], fields: list[str], delimiter: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def split_field(value: str) -> list[str]:
    return [part.strip() for part in (value or "").split(";") if part.strip()]


def ena_url(study: str) -> str:
    return (
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + study
        + "&result=read_run&fields=" + ENA_FIELDS + "&format=tsv&download=false"
    )


def selected_rows(study: str) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    source = ROOT / "provenance" / study / "source_snapshots" / "ena_read_run.tsv"
    with source.open(encoding="utf-8-sig", newline="") as handle:
        all_rows = list(csv.DictReader(handle, delimiter="\t"))
    selected, excluded = [], []
    for row in all_rows:
        keep = (
            row.get("library_strategy") == "RNA-Seq"
            and row.get("scientific_name", "").startswith("Schistosoma mansoni")
            and bool(row.get("fastq_ftp"))
        )
        (selected if keep else excluded).append(row)
    return selected, excluded


def parse_design(study: str, row: dict[str, str]) -> dict[str, str]:
    name = row.get("library_name") or row.get("sample_title") or row["run_accession"]
    base = {
        "source_sample_name": name,
        "biosample": row.get("sample_accession", ""),
        "stage": "adult",
        "tissue": "whole_organism",
        "sex": "mixed",
        "infection_mode": "not_declared",
        "time_hours": "not_applicable",
        "batch": "not_declared",
    }
    if study == "PRJNA312093":
        lower = name.lower()
        sex = "female" if "female" in lower else "male"
        rep = re.search(r"replicate\s+(\d+)", lower).group(1)
        if "cercaria" in lower:
            stage = "cercaria"
        elif "stage 1" in lower:
            stage = "schistosomulum_stage_1"
        elif "stage 2" in lower:
            stage = "schistosomulum_stage_2"
        elif "stage 3" in lower:
            stage = "schistosomulum_stage_3"
        elif "adulte" in lower or "adult" in lower:
            stage = "adult"
        else:
            raise ValueError(f"{study}: estágio não reconhecido em {name!r}")
        return {**base, "sample_id": f"{sex}_{stage}_r{rep}", "condition": f"{stage}_{sex}", "stage": stage, "sex": sex, "treatment": "untreated", "replicate": rep, "technical_unit": "run_1"}
    if study == "PRJNA494090":
        match = re.fullmatch(r"P(\d+)_(control|PZQ-selected)(\d+)", row.get("sample_title", ""), re.I)
        if not match:
            raise ValueError(f"{study}: nome não reconhecido em {row.get('sample_title')!r}")
        passage, group, rep = match.groups()
        selected = group.lower() != "control"
        condition = f"p{passage}_{'pzq_selected' if selected else 'control'}"
        return {**base, "sample_id": f"{condition}_r{rep}", "condition": condition, "stage": "not_declared", "treatment": "pzq_selected_line" if selected else "control_line", "batch": f"passage_{passage}", "replicate": rep, "technical_unit": "run_1"}
    if study == "PRJNA704646":
        match = re.fullmatch(r"SmLE-PZQ-(ER|ES)-RNA-(adu|juv)-(m|f)\.(\d+)", name)
        if not match:
            raise ValueError(f"{study}: nome não reconhecido em {name!r}")
        allele, stage_code, sex_code, library_rep = match.groups()
        stage = {"adu": "adult", "juv": "juvenile"}[stage_code]
        sex = {"m": "male", "f": "female"}[sex_code]
        allele_name = {"ER": "resistant_enriched", "ES": "sensitive_enriched"}[allele]
        condition = f"{allele_name}_{stage}_{sex}"
        return {**base, "sample_id": condition, "condition": condition, "stage": stage, "sex": sex, "treatment": allele_name, "replicate": "1", "technical_unit": f"library_rep_{library_rep}"}
    if study == "PRJNA1150591":
        match = re.fullmatch(r"(4wk|7wk)_(MCLZ|CTRL)_([A-E])", name)
        if not match:
            raise ValueError(f"{study}: nome não reconhecido em {name!r}")
        week, group, letter = match.groups()
        rep = str(ord(letter) - ord("A") + 1)
        stage = "juvenile_4wk" if week == "4wk" else "adult_7wk"
        treatment = "meclonazepam" if group == "MCLZ" else "control"
        condition = f"{stage}_{'mclz' if group == 'MCLZ' else 'control'}"
        return {**base, "sample_id": f"{condition}_r{rep}", "condition": condition, "stage": stage, "treatment": treatment, "batch": week, "replicate": rep, "technical_unit": "run_1"}
    if study == "PRJNA1177852":
        match = re.fullmatch(r"MCLZ_(15min|\d+hr)", name)
        if not match:
            raise ValueError(f"{study}: nome não reconhecido em {name!r}")
        label = match.group(1)
        hours = "0.25" if label == "15min" else label.removesuffix("hr")
        condition = "mclz_time_" + hours.replace(".", "p") + "h"
        return {**base, "sample_id": condition, "condition": condition, "treatment": "meclonazepam", "time_hours": hours, "replicate": "1", "technical_unit": "run_1"}
    if study == "PRJNA551487":
        match = re.fullmatch(r"(FPL|DMSO)_(\d+)", name)
        if not match:
            raise ValueError(f"{study}: nome não reconhecido em {name!r}")
        group, rep = match.groups()
        condition = "adult_fpl64176" if group == "FPL" else "adult_dmso"
        return {**base, "sample_id": f"{condition}_r{rep}", "condition": condition, "treatment": "fpl_64176" if group == "FPL" else "dmso", "replicate": rep, "technical_unit": "run_1"}
    if study == "PRJNA428470":
        match = re.fullmatch(r"(5-AzaC|Control)(\d+)", name)
        if not match:
            raise ValueError(f"{study}: nome não reconhecido em {name!r}")
        group, rep = match.groups()
        condition = "female_adult_5azac" if group == "5-AzaC" else "female_adult_control"
        return {**base, "sample_id": f"{condition}_r{rep}", "condition": condition, "sex": "female", "treatment": "5_azacitidine" if group == "5-AzaC" else "control", "replicate": rep, "technical_unit": "run_1"}
    raise ValueError(f"Estudo sem parser: {study}")


def de_spec(study: str) -> dict[str, object] | None:
    contrast_levels = {
        "PRJNA312093": [(f"{stage}_female", f"{stage}_male") for stage in ("cercaria", "schistosomulum_stage_1", "schistosomulum_stage_2", "schistosomulum_stage_3", "adult")],
        "PRJNA494090": [(f"p{p}_pzq_selected", f"p{p}_control") for p in (6, 7, 8)],
        "PRJNA1150591": [("adult_7wk_mclz", "adult_7wk_control")],
        "PRJNA551487": [("adult_fpl64176", "adult_dmso")],
        "PRJNA428470": [("female_adult_5azac", "female_adult_control")],
    }
    if study not in contrast_levels:
        return None
    contrasts = []
    for numerator, denominator in contrast_levels[study]:
        contrasts.append({
            "id": f"condition__{numerator}_vs_{denominator}",
            "factor": "condition", "numerator": numerator, "denominator": denominator,
            "description": f"{numerator} versus {denominator}",
            "direction": f"{numerator}/{denominator}",
        })
    return {
        "schema_version": "1.0", "scope": study, "correction": "raw", "provider": "deseq2", "test": "wald",
        "design": {"variable": "condition", "covariates": [], "formula": "~ condition"},
        "filter": {"method": "total_count", "operator": ">=", "threshold": 10},
        "parameters": {"alpha": 0.05, "lfc_threshold": 1.0, "min_replicates": 2, "non_integer_counts": "round"},
        "analysis_id": f"smansoni_{study.lower()}_primary",
        "target_dir": f"results/rnaseq/differential_expression/{study}",
        "contrasts": contrasts,
    }


def blocked_spec(study: str) -> dict[str, object]:
    if study == "PRJNA704646":
        code = "NO_BIOLOGICAL_REPLICATION_AFTER_TECHNICAL_CONSOLIDATION"
        reason = "Três bibliotecas derivam do mesmo pool de RNA em cada combinação; após a consolidação há n=1 biológico por grupo."
    else:
        code = "UNREPLICATED_TIMECOURSE_WITHOUT_PARALLEL_CONTROL"
        reason = "Há uma biblioteca por tempo e nenhum controle paralelo replicado; a série é apenas descritiva."
    return {
        "schema_version": "1.0", "scope": study, "status": "BLOCKED", "reason_code": code,
        "reason": reason,
        "allowed_outputs": ["FASTQC", "TRIM_GALORE", "SALMON_QUANT", "TXIMPORT_COUNTS", "TPM"],
        "forbidden_outputs": ["DESEQ2_INFERENCE"],
    }


def make_study(study: str) -> dict[str, object]:
    spec = STUDIES[study]
    raw_rows, excluded = selected_rows(study)
    if len(raw_rows) != spec["expected_runs"]:
        raise ValueError(f"{study}: esperados {spec['expected_runs']} runs, encontrados {len(raw_rows)}")
    parsed = [(row, parse_design(study, row)) for row in raw_rows]
    run_rows: list[dict[str, str]] = []
    metadata_rows: list[dict[str, str]] = []
    grouped: defaultdict[str, list[tuple[dict[str, str], dict[str, str]]]] = defaultdict(list)
    for row, design in parsed:
        urls = split_field(row.get("fastq_ftp", ""))
        md5s = split_field(row.get("fastq_md5", ""))
        sizes = split_field(row.get("fastq_bytes", ""))
        expected_files = 1 if spec["layout"] == "SINGLE" else 2
        if not (len(urls) == len(md5s) == len(sizes) == expected_files):
            raise ValueError(f"{study}/{row['run_accession']}: conjunto FASTQ incompleto")
        urls = [("https://" + item) if item.startswith("ftp.") else item for item in urls]
        run_rows.append({
            "dataset": study, "sample_id": design["sample_id"], "run_accession": row["run_accession"],
            "biosample": row.get("sample_accession", ""), "secondary_biosample": row.get("secondary_sample_accession", ""),
            "experiment_accession": row.get("experiment_accession", ""), "technical_unit": design["technical_unit"],
            "library_name": row.get("library_name", "") or row.get("sample_title", ""), "library_layout": spec["layout"],
            "fastq_1_url": urls[0], "fastq_2_url": urls[1] if expected_files == 2 else "",
            "fastq_1_md5": md5s[0], "fastq_2_md5": md5s[1] if expected_files == 2 else "",
            "fastq_1_bytes": sizes[0], "fastq_2_bytes": sizes[1] if expected_files == 2 else "",
            "read_count": row.get("read_count", ""), "base_count": row.get("base_count", ""),
            "ena_first_public": row.get("first_public", ""), "ena_query_url": ena_url(study),
        })
        metadata_rows.append({
            "dataset": study, "sample_id": design["sample_id"], "file_prefix": design["sample_id"],
            "run_accession": row["run_accession"], "condition": design["condition"], "stage": design["stage"],
            "tissue": design["tissue"], "sex": design["sex"], "infection_mode": design["infection_mode"],
            "treatment": design["treatment"], "time_hours": design["time_hours"], "batch": design["batch"],
            "replicate": design["replicate"], "biosample": row.get("sample_accession", ""),
            "library_name": row.get("library_name", "") or row.get("sample_title", ""), "technical_unit": design["technical_unit"],
            "library_layout": spec["layout"].lower(),
            "fragment_length_mean": spec.get("fragment_mean", ""), "fragment_length_sd": spec.get("fragment_sd", ""),
        })
        grouped[design["sample_id"]].append((row, design))

    sample_rows = []
    for sample_id, members in sorted(grouped.items()):
        designs = [item[1] for item in members]
        first = designs[0]
        invariant = ("condition", "stage", "tissue", "sex", "infection_mode", "treatment", "time_hours", "batch", "replicate")
        for field in invariant:
            if len({d[field] for d in designs}) != 1:
                raise ValueError(f"{study}/{sample_id}: {field} inconsistente")
        sample_rows.append({
            "dataset": study, "sample_id": sample_id, "source_sample_name": first["source_sample_name"],
            "biosample": ";".join(sorted({item[0].get("sample_accession", "") for item in members})),
            **{field: first[field] for field in invariant}, "technical_runs": str(len(members)),
        })
    if len(sample_rows) != spec["expected_samples"]:
        raise ValueError(f"{study}: esperadas {spec['expected_samples']} amostras, encontradas {len(sample_rows)}")

    run_rows.sort(key=lambda r: (r["sample_id"], r["run_accession"]))
    metadata_rows.sort(key=lambda r: (r["sample_id"], r["run_accession"]))
    write_rows(ROOT / "metadata" / study / "runs.tsv", run_rows, RUN_FIELDS, "\t")
    write_rows(ROOT / "metadata" / study / "samples.tsv", sample_rows, SAMPLE_FIELDS, "\t")
    write_rows(ROOT / "metadata" / study / "metadata.csv", metadata_rows, METADATA_FIELDS, ",")
    if excluded:
        excluded_rows = [
            {"run_accession": row.get("run_accession", ""), "library_strategy": row.get("library_strategy", ""),
             "scientific_name": row.get("scientific_name", ""), "reason": "NOT_TARGET_RNASEQ_S_MANSONI_OR_FASTQ_UNAVAILABLE"}
            for row in excluded
        ]
        write_rows(ROOT / "metadata" / study / "excluded_ena_runs.tsv", excluded_rows,
                   ["run_accession", "library_strategy", "scientific_name", "reason"], "\t")

    config_dir = ROOT / "config" / study
    config_dir.mkdir(parents=True, exist_ok=True)
    (config_dir / "pipeline_config.sh").write_text(
        '#!/usr/bin/env bash\n: "${HF_HELIXFORGE_ROOT:?set HF_HELIXFORGE_ROOT}"\n: "${HF_PACKAGE_ROOT:?set HF_PACKAGE_ROOT}"\n'
        f'export PROJECT_DIR="${{HF_HELIXFORGE_ROOT}}/pipelines/rnaseq"\nexport USER_SETTINGS_FILE="${{HF_PACKAGE_ROOT}}/config/{study}/user_settings.sh"\n'
        '# shellcheck disable=SC1090\nsource "${HF_HELIXFORGE_ROOT}/pipelines/rnaseq/config/pipeline_config.sh"\n',
        encoding="utf-8",
    )
    settings = f'''#!/usr/bin/env bash
: "${{HF_PACKAGE_ROOT:?source config/server.env before running}}"
: "${{HF_DATA_ROOT:?set HF_DATA_ROOT}}"
: "${{HF_REFERENCE_ROOT:?set HF_REFERENCE_ROOT}}"
: "${{HF_RESULTS_ROOT:?set HF_RESULTS_ROOT}}"

export PIPELINE_NAME="smansoni_{study.lower()}_wbps19"
export ORGANISM_NAME="Schistosoma_mansoni"
export REFERENCE_ID="Schistosoma_mansoni_SM_V10_WBPS19"
export PIPELINE_PROJECTS="{study}"
export SCRATCH_ROOT="${{HF_DATA_ROOT}}"
export METADATA_FINAL="${{HF_PACKAGE_ROOT}}/metadata/{study}/metadata.csv"
export METADATA_FINAL_NEW="${{METADATA_FINAL}}"
export REF_GENOME_FA="${{HF_REFERENCE_ROOT}}/SM_V10_WBPS19/schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa"
export REF_TRANSCRIPTS_FA="${{HF_REFERENCE_ROOT}}/SM_V10_WBPS19/schistosoma_mansoni.PRJEA36577.WBPS19.mRNA_transcripts.fa"
export REF_GTF="${{HF_REFERENCE_ROOT}}/SM_V10_WBPS19/schistosoma_mansoni.PRJEA36577.WBPS19.canonical_geneset.gtf"
export REF_GFF3="${{HF_REFERENCE_ROOT}}/SM_V10_WBPS19/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3"
export SALMON_INDEX_DIR="${{HF_REFERENCE_ROOT}}/SM_V10_WBPS19/salmon_index_k31"
export QUANT_DIR="${{HF_RESULTS_ROOT}}/{study}/compatibility/quants"
export QUANT_METHOD="salmon"
export SALMON_KMER_SIZE=31
export TRIM_QUALITY=20
export TRIM_LENGTH=20
export RUN_BATCH_CORRECTION=0
export RUN_GENE_REPORT=0
export FASTQ_LAYOUT="{spec['layout'].lower()}"
export HF_RNASEQ_RUN_MODE="{spec['mode']}"
'''
    (config_dir / "user_settings.sh").write_text(settings, encoding="utf-8")
    if (document := de_spec(study)) is not None:
        (config_dir / "de_spec.json").write_text(json.dumps(document, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    else:
        (config_dir / "de_spec.blocked.json").write_text(json.dumps(blocked_spec(study), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    total_bytes = sum(int(item) for row in run_rows for item in (row["fastq_1_bytes"], row["fastq_2_bytes"]) if item)
    readme = f"""# {study}

## Decisão científica

{spec['title']}. Fonte primária: {spec['paper']}.

Foram selecionados **{len(run_rows)} runs RNA-seq** de *Schistosoma mansoni*, organizados em **{len(sample_rows)} amostras biológicas**. Layout: `{spec['layout']}`. Modo HelixForge: `{spec['mode']}`.

## Limitação congelada

{spec['limitations']}
"""
    if study == "PRJNA312093":
        readme += f"\nPrior single-end: média `{spec['fragment_mean']}` e desvio `{spec['fragment_sd']}`. Fonte: {spec['fragment_source']}. Esses valores não foram relatados no artigo e devem passar por análise de sensibilidade antes da incorporação ao atlas.\n"
    readme += "\nOs FASTQs brutos não são versionados. `runs.tsv` fixa URLs, MD5, tamanhos e a relação run → amostra; o snapshot ENA imutável está em `provenance/%s/source_snapshots/ena_read_run.tsv`.\n" % study
    (config_dir / "README.md").write_text(readme, encoding="utf-8")

    source = ROOT / "provenance" / study / "source_snapshots" / "ena_read_run.tsv"
    manifest = {
        "schema_version": "1.0", "study": study, "prepared_on": PREPARED_ON,
        "reference": "Schistosoma mansoni SM_V10 / WormBase ParaSite WBPS19",
        "ena_query_url": ena_url(study),
        "source_snapshot": {"path": str(source.relative_to(ROOT)).replace("\\", "/"), "sha256": sha256(source), "bytes": source.stat().st_size},
        "selection": {"library_strategy": "RNA-Seq", "organism_prefix": "Schistosoma mansoni", "fastq_required": True,
                      "selected_runs": len(run_rows), "excluded_rows": len(excluded)},
        "logical_samples": len(sample_rows), "library_layout": spec["layout"],
        "raw_fastq_bytes": total_bytes, "raw_fastq_decimal_gb": round(total_bytes / 1_000_000_000, 3),
        "run_mode": spec["mode"], "differential_expression": "ENABLED" if document is not None else "BLOCKED",
    }
    if study == "PRJNA312093":
        manifest["single_end_fragment_prior"] = {"mean": 200, "sd": 80, "source": spec["fragment_source"], "status": "PROVISIONAL_BIOLOGICAL_GATE"}
    prov_path = ROOT / "provenance" / study / "input_manifest.json"
    prov_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    manifests = [make_study(study) for study in STUDIES]
    pending_rows = []
    for manifest in manifests:
        pending_rows.append({
            "study": manifest["study"], "priority": "P1", "preparation_status": "READY_FOR_SERVER_PREFLIGHT",
            "run_mode": manifest["run_mode"], "library_layout": manifest["library_layout"],
            "runs": str(manifest["selection"]["selected_runs"]), "biological_samples": str(manifest["logical_samples"]),
            "raw_fastq_bytes": str(manifest["raw_fastq_bytes"]), "differential_expression": manifest["differential_expression"],
            "metadata_path": f"metadata/{manifest['study']}/metadata.csv", "provenance_path": f"provenance/{manifest['study']}/input_manifest.json",
        })
    write_rows(ROOT / "metadata" / "p1_pending_manifest.tsv", pending_rows,
               ["study", "priority", "preparation_status", "run_mode", "library_layout", "runs", "biological_samples", "raw_fastq_bytes", "differential_expression", "metadata_path", "provenance_path"], "\t")
    global_provenance = {
        "schema_version": "1.0", "prepared_on": PREPARED_ON,
        "scope": "P1_INPUT_PREPARATION_ONLY", "atlas_integration": "NOT_PERFORMED",
        "execution_plan": "LOCAL_ONLY_NOT_VERSIONED_IN_THIS_REPOSITORY",
        "studies": manifests,
        "totals": {"studies": len(manifests), "runs": sum(m["selection"]["selected_runs"] for m in manifests),
                   "biological_samples": sum(m["logical_samples"] for m in manifests),
                   "raw_fastq_bytes": sum(m["raw_fastq_bytes"] for m in manifests)},
    }
    (ROOT / "provenance" / "p1_pending_preparation.json").write_text(json.dumps(global_provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(global_provenance["totals"], indent=2))


if __name__ == "__main__":
    main()
