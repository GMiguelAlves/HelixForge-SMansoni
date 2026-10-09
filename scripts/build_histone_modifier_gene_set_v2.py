#!/usr/bin/env python3
"""Build the versioned SM_V10 histone-modifier catalogue.

The historical list is treated as a user-curated biological classification.
Identifier resolution against WBPS19/SM_V10 is kept separate from biological
curation so that an ambiguous legacy ID is never mislabeled as weak evidence.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "config" / "PRJEB32839" / "candidate_genes.txt"
ANNOTATIONS = ROOT / "config" / "reference" / "gene_annotations_wbps19.tsv"
OUTDIR = ROOT / "config" / "candidates"
TSV_OUT = OUTDIR / "histone_modifier_gene_set_v2.tsv"
DOC_OUT = ROOT / "docs" / "conjunto_modificadores_histonas_v2.md"
PROV_OUT = ROOT / "provenance" / "histone_modifier_gene_set_v2.json"
VALIDATION_TABLE = ROOT / "provenance" / "histone_modifier_review_validation" / "review_add_validation.tsv"
VALIDATION_MANIFEST = ROOT / "provenance" / "histone_modifier_review_validation" / "validation_manifest.json"
ID_RESOLUTION_TABLE = ROOT / "provenance" / "histone_modifier_id_resolution" / "historical_id_resolution.tsv"
ID_RESOLUTION_MANIFEST = ROOT / "provenance" / "histone_modifier_id_resolution" / "validation_manifest.json"

CATALOGUE_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6125318/"
BROMO_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC10726355/"


ID_MAP = {
    "Smp_064830": "Smp_341020",
    "Smp_081570": "Smp_246840",
    "Smp_097090": "Smp_341630",
    "Smp_124950": "Smp_342100",
    "Smp_128500": "Smp_342360",
    "Smp_131300": "Smp_016750",
    "Smp_137060": "Smp_307060",
    "Smp_154860": "Smp_344680",
    "Smp_174840": "Smp_346320",
    "Smp_196270": "Smp_315890",
    "Smp_211290": "Smp_337860",
}

VALIDATED_ID_MAP = {
    "Smp_137240": ("Smp_316180", "RESOLVED_SPLIT_FUNCTIONAL_SUCCESSOR", "Locus histórico fragmentado; Smp_316180 preserva o segmento JmjC funcional."),
    "Smp_147870": ("Smp_213920", "RESOLVED_SPLIT_FUNCTIONAL_SUCCESSOR", "Locus histórico dividido; Smp_213920 preserva o segmento JmjC, enquanto Smp_213930 corresponde ao segmento PPM."),
    "Smp_180990": ("Smp_241580", "RESOLVED_BY_SEQUENCE_EXACT", "Sequência histórica mapeada com 100% de identidade e 100% de cobertura."),
    "Smp_161010": ("Smp_246410", "RESOLVED_BY_SEQUENCE_HIGH_CONFIDENCE", "Melhor correspondência SM_V10 com 98,637% de identidade e 99% de cobertura."),
    "Smp_163470": ("Smp_246920", "RESOLVED_BY_SEQUENCE_EXACT", "Sequência histórica mapeada com 100% de identidade e 100% de cobertura."),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_annotations() -> dict[str, dict[str, str]]:
    with ANNOTATIONS.open(encoding="utf-8", newline="") as handle:
        return {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def read_historical_groups() -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {}
    for raw in SOURCE.read_text(encoding="utf-8").splitlines():
        if not raw.strip():
            continue
        label, genes = raw.split(":", 1)
        groups[label.strip()] = list(dict.fromkeys(re.findall(r"Smp_\d+", genes)))
    return groups


def resolve(gene: str, annotations: dict[str, dict[str, str]]) -> tuple[str, str, str]:
    if gene in annotations:
        return gene, "EXACT_SM_V10", ""
    if gene in ID_MAP:
        return ID_MAP[gene], "UNIQUE_PREVIOUS_STABLE_ID", ""
    if gene in VALIDATED_ID_MAP:
        return VALIDATED_ID_MAP[gene]
    return "", "NOT_FOUND_IN_WBPS19", "Revisar mapeamento de identificador."


def main() -> None:
    annotations = read_annotations()
    historical = read_historical_groups()
    rows: list[dict[str, str]] = []

    def add(
        input_gene_id: str,
        broad_class: str,
        functional_group: str,
        scope_tier: str,
        curation_status: str,
        evidence_source: str,
        source_url: str = "",
        notes: str = "",
    ) -> None:
        current, mapping_status, mapping_note = resolve(input_gene_id, annotations)
        names = []
        if mapping_status in {"EXACT_SM_V10", "UNIQUE_PREVIOUS_STABLE_ID"} or mapping_status.startswith("RESOLVED_"):
            names.append(annotations[current]["functional_name"])
        rows.append(
            {
                "input_gene_id": input_gene_id,
                "sm_v10_gene_id": current,
                "broad_class": broad_class,
                "functional_group": functional_group,
                "scope_tier": scope_tier,
                "curation_status": curation_status,
                "mapping_status": mapping_status,
                "functional_name_wbps19": " | ".join(dict.fromkeys(names)),
                "evidence_source": evidence_source,
                "source_url": source_url,
                "notes": " ".join(x for x in (notes, mapping_note) if x).strip(),
            }
        )

    historical_spec = {
        "Grupo HP1": ("HP1", "HETEROCHROMATIN_HP1", "CONTEXT"),
        "Writers": ("WRITER", "HISTONE_METHYLATION_WRITER_CURATED", "CORE_CURATED"),
        "Readers": ("READER", "HISTONE_METHYLATION_READER_CURATED", "CORE_CURATED"),
        "Erasers": ("ERASER", "HISTONE_METHYLATION_ERASER_CURATED", "CORE_CURATED"),
    }
    for label in ("Grupo HP1", "Writers", "Readers", "Erasers"):
        broad, group, tier = historical_spec[label]
        for gene in historical[label]:
            note = (
                "Grupo mantido à parte: HP1 é leitor/organizador de heterocromatina, "
                "mas não é contado como modificador enzimático no núcleo."
                if label == "Grupo HP1"
                else "Classificação biológica confirmada pela curadoria do projeto."
            )
            add(
                gene,
                broad,
                group,
                tier,
                "USER_CURATED_CONFIRMED",
                "USER_CURATED_HISTORICAL_LIST",
                "",
                note,
            )

    explicit_groups = {
        ("WRITER", "HISTONE_ACETYLTRANSFERASE", "CORE_ADD"): [
            "Smp_053140", "Smp_070190", "Smp_105910", "Smp_127010",
            "Smp_131320", "Smp_178700", "Smp_194520", "Smp_198670",
        ],
        ("ERASER", "HISTONE_DEACETYLASE", "CORE_ADD"): [
            "Smp_005210", "Smp_069380", "Smp_091990", "Smp_093280",
            "Smp_138770", "Smp_191310",
        ],
        ("ERASER", "SIRTUIN_DEACETYLASE", "CORE_ADD"): ["Smp_210340", "Smp_340470"],
        ("WRITER", "HISTONE_LYSINE_METHYLTRANSFERASE", "CORE_ADD"): ["Smp_056930", "Smp_247270"],
    }
    for (broad, group, tier), genes in explicit_groups.items():
        for gene in genes:
            add(
                gene,
                broad,
                group,
                tier,
                "SM_V10_EXPLICIT_ANNOTATION",
                "WORMBASE_PARASITE_WBPS19_SM_V10",
                notes="Ausente da lista histórica e sustentado por nome funcional explícito em WBPS19/SM_V10.",
            )

    add(
        "Smp_084140",
        "ERASER",
        "SIRTUIN_DEACETYLASE",
        "CORE_ADD",
        "VALIDATED_SIR2_DOMAIN_AND_SIRT2_ORTHOLOGY",
        "WBPS19_UNIPROT_PFAM_BLASTP_VALIDATION",
        "https://rest.uniprot.org/uniprotkb/A0A5K4EHC2",
        "PF02146/SIR2 completo; melhor hit humano SIRT2, 46,1% identidade, 89% de cobertura, E=5,46e-97.",
    )
    add(
        "Smp_127180",
        "READER",
        "TAF3_PHD_H3K4ME3_READER",
        "CORE_ADD",
        "VALIDATED_TAF3_PHD_READER_ORTHOLOGY",
        "WBPS19_ORTHOLOGY_PFAM_BLASTP_AND_TAF3_LITERATURE",
        "https://pubmed.ncbi.nlm.nih.gov/18682226/",
        "PHD PF00628 e ortologia TAF3; o PHD de TAF3 reconhece H3K4me3.",
    )

    review_groups = {
        ("READER", "PHD_READER_CANDIDATE"): ["Smp_185180", "Smp_316920"],
    }
    for (broad, group), genes in review_groups.items():
        for gene in genes:
            add(
                gene,
                broad,
                group,
                "REVIEW_ADD",
                "DOMAIN_ORTHOLOGY_REVIEW_REQUIRED",
                "WBPS19_PFAM_BLASTP_VALIDATION",
                notes=(
                    "Domínio PHD confirmado, mas a especificidade por histona não foi estabelecida. "
                    "Smp_316920 é um modelo curto com PHDs em tandem; revisar o modelo gênico."
                    if gene == "Smp_316920"
                    else "Domínio PHD confirmado; sem ortologia funcional suficiente para promover ao núcleo."
                ),
            )

    excluded_review = [
        (
            "Smp_024380", "WRITER", "NON_HISTONE_RMT2_GAMT_LIKE", "EXCLUDED_NON_HISTONE_METHYLTRANSFERASE",
            "O domínio RMT2 também ocorre em GAMT; o melhor hit humano é GAMT e não há evidência de metilação de histona.",
        ),
        (
            "Smp_059920", "WRITER", "EZH_ANNOTATION_CONFLICT", "EXCLUDED_NO_SET_DOMAIN",
            "Proteína de 210 aa sem domínio SET; ortólogo de S. haematobium foi anotado como proteína com domínio B3.",
        ),
        (
            "Smp_064250", "READER", "PHD_TRANSCRIPTION_FACTOR_CONTEXT", "EXCLUDED_NO_HISTONE_READER_ORTHOLOGY",
            "PHD confirmado, mas a ortologia dominante é TCF20; não há suporte para leitura direta de marca de histona.",
        ),
        (
            "Smp_128520", "ERASER", "JMJC_ANNOTATION_CONFLICT", "EXCLUDED_NO_CANONICAL_JMJC_DOMAIN",
            "Sem hit PF02373, sem motivo catalítico HxD e sem ortólogo humano detectável no BLASTP.",
        ),
        (
            "Smp_128530", "ERASER", "JMJC_ANNOTATION_CONFLICT", "EXCLUDED_NO_CANONICAL_JMJC_DOMAIN",
            "Sem hit PF02373, sem motivo catalítico HxD e sem ortólogo humano detectável no BLASTP.",
        ),
        (
            "Smp_204920", "READER", "TUDOR_ANNOTATION_CONFLICT", "EXCLUDED_NO_CANONICAL_TUDOR_DOMAIN",
            "Sem hit PF00567 e sem ortologia funcional que sustente leitura de histona.",
        ),
        (
            "Smp_329500", "READER", "TUDOR_ANNOTATION_CONFLICT", "EXCLUDED_SHORT_UNSUPPORTED_MODEL",
            "Modelo de 70 aa sem hit PF00567 e sem ortologia informativa; inadequado para o núcleo.",
        ),
    ]
    for gene, broad, group, status, note in excluded_review:
        add(
            gene,
            broad,
            group,
            "EXCLUDED_REVIEW",
            status,
            "WBPS19_UNIPROT_PFAM_BLASTP_VALIDATION",
            notes=note,
        )

    bromodomain_genes = [
        "Smp_070190", "Smp_267320", "Smp_147950", "Smp_105910", "Smp_122370",
        "Smp_127010", "Smp_131320", "Smp_337800", "Smp_168890", "Smp_179290",
        "Smp_246920", "Smp_170760", "Smp_170390", "Smp_213110", "Smp_345670",
        "Smp_128120", "Smp_158050", "Smp_164080", "Smp_082660", "Smp_083790",
        "Smp_096240", "Smp_094380",
    ]
    for gene in bromodomain_genes:
        add(
            gene,
            "READER",
            "BROMODOMAIN_ACETYL_READER",
            "CORE_ADD",
            "PUBLISHED_BROMODOMAIN_CATALOGUE",
            "BROMODOMAIN_CATALOGUE_2023",
            BROMO_URL,
            "Catálogo publicado de proteínas com bromodomínio; genes multifuncionais mantêm associações adicionais.",
        )

    contextual_groups = {
        ("COFACTOR", "POLYCOMB_TRITHORAX_COREPRESSOR_CONTEXT"): [
            "Smp_014400", "Smp_035290", "Smp_064650", "Smp_068950", "Smp_085480",
            "Smp_089230", "Smp_142170", "Smp_212260", "Smp_330000",
        ],
        ("COFACTOR", "TRANSCRIPTIONAL_PHD_CONTEXT"): ["Smp_064250"],
        ("REMODELER", "ATP_DEPENDENT_CHROMATIN_REMODELING"): [
            "Smp_103530", "Smp_148010", "Smp_152650", "Smp_154640",
            "Smp_158050", "Smp_163620", "Smp_243660", "Smp_267320",
        ],
        ("CHAPERONE", "HISTONE_CHAPERONE_NUCLEOSOME_ASSEMBLY"): [
            "Smp_153610", "Smp_169140", "Smp_320670", "Smp_333660",
        ],
    }
    for (broad, group), genes in contextual_groups.items():
        for gene in genes:
            note = "Maquinaria cromatínica contextual; não contar como modificador direto."
            if gene in {"Smp_169140", "Smp_320670"}:
                note += " Revisar relação entre loci antes de agregação gene-cêntrica."
            add(
                gene,
                broad,
                group,
                "CONTEXT",
                "SM_V10_CONTEXTUAL_ANNOTATION",
                "WORMBASE_PARASITE_WBPS19_SM_V10",
                notes=note,
            )

    rank = {"CORE_CURATED": 0, "CORE_ADD": 1, "REVIEW_ADD": 2, "CONTEXT": 3, "EXCLUDED_REVIEW": 4}
    rows.sort(key=lambda r: (rank[r["scope_tier"]], r["broad_class"], r["functional_group"], r["input_gene_id"]))
    for number, row in enumerate(rows, start=1):
        row["record_id"] = f"HMEV2-{number:04d}"

    fields = [
        "record_id", "input_gene_id", "sm_v10_gene_id", "broad_class", "functional_group",
        "scope_tier", "curation_status", "mapping_status", "functional_name_wbps19",
        "evidence_source", "source_url", "notes",
    ]
    OUTDIR.mkdir(parents=True, exist_ok=True)
    with TSV_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    unique_by_tier = {
        tier: len({r["input_gene_id"] for r in rows if r["scope_tier"] == tier})
        for tier in rank
    }
    unique_all = len({r["input_gene_id"] for r in rows})
    class_counts = Counter((r["scope_tier"], r["broad_class"]) for r in rows)
    mapping_counts = Counter(
        r["mapping_status"] for r in rows if r["evidence_source"] == "USER_CURATED_HISTORICAL_LIST"
    )

    doc = f"""# Conjunto de modificadores de histonas v2 (SM_V10/WBPS19)

## Decisão de curadoria

A lista histórica é uma curadoria biológica do projeto. Os membros de `Writers`,
`Readers` e `Erasers` são tratados como modificadores de histonas confirmados pela
curadoria da autora. O grupo HP1 é preservado, mas separado como leitor/organizador
de heterocromatina e não entra na contagem de modificadores enzimáticos.

O estado do identificador é uma questão diferente da classificação biológica. Os cinco
IDs históricos que antes estavam ambíguos ou ausentes foram resolvidos por comparação
de sequência e, nos loci fragmentados, pela preservação do segmento funcional esperado.

## Camadas da lista

| Camada | Critério | Genes únicos | Uso no atlas |
|---|---|---:|---|
| `CORE_CURATED` | lista histórica curada pelo projeto, excluindo HP1 da contagem enzimática | {unique_by_tier['CORE_CURATED']} | núcleo, após resolver IDs pendentes |
| `CORE_ADD` | anotação explícita, catálogo publicado ou validação de domínio+ortologia | {unique_by_tier['CORE_ADD']} | expansão principal |
| `REVIEW_ADD` | domínio confirmado, mas função sobre histona ou modelo gênico ainda inconclusivo | {unique_by_tier['REVIEW_ADD']} | painel exploratório, fora do núcleo |
| `CONTEXT` | HP1, cofatores, remodeladores e chaperonas | {unique_by_tier['CONTEXT']} | contexto mecanístico, fora da contagem de modificadores diretos |
| `EXCLUDED_REVIEW` | hipótese inicial não sustentada como modificador/reader direto | {unique_by_tier['EXCLUDED_REVIEW']} | documentado, mas não analisado como candidato |

Há **{unique_all} genes de entrada únicos** em todas as camadas. A tabela é
normalizada por associação funcional: um gene pode aparecer em mais de uma linha,
como uma acetiltransferase que também contém bromodomínio.

## Classes importantes acrescentadas

- **Writers de acetilação:** HATs explícitas em SM_V10, ausentes da lista histórica.
- **Erasers de acetilação:** HDACs clássicas e sirtuínas; `Smp_084140` foi
  promovido após confirmação do domínio SIR2 e ortologia com SIRT2.
- **Readers de acetilação:** 22 proteínas do catálogo publicado de bromodomínios.
- **Reader de H3K4me3:** `Smp_127180` foi promovido como ortólogo de TAF3 com
  PHD; a função de leitura de H3K4me3 do PHD de TAF3 é experimentalmente estabelecida.
- **PHDs ainda exploratórios:** `Smp_185180` e `Smp_316920` têm domínio PHD,
  mas ainda não possuem suporte suficiente para especificidade de marca de histona.
- **Hipóteses rejeitadas para o painel direto:** `Smp_024380`, `Smp_059920`,
  `Smp_064250`, `Smp_128520`, `Smp_128530`, `Smp_204920` e `Smp_329500` foram
  mantidos como registros de auditoria em `EXCLUDED_REVIEW`.
- **Contexto cromatínico:** Polycomb/Trithorax, CoREST/SAP18/RBBP4, remodeladores
  ATP-dependentes e chaperonas de histonas são mantidos em grupos funcionais próprios.

## Contagens por camada e classe

| Camada | Classe | Associações gene-classe |
|---|---|---:|
"""
    for (tier, broad), count in sorted(class_counts.items(), key=lambda x: (rank[x[0][0]], x[0][1])):
        doc += f"| `{tier}` | `{broad}` | {count} |\n"

    doc += f"""

## Estado dos IDs da curadoria histórica

- correspondência exata SM_V10: {mapping_counts['EXACT_SM_V10']} associações;
- correspondência única por ID anterior: {mapping_counts['UNIQUE_PREVIOUS_STABLE_ID']} associações;
- correspondência exata por sequência: {mapping_counts['RESOLVED_BY_SEQUENCE_EXACT']} associações;
- correspondência de alta confiança por sequência: {mapping_counts['RESOLVED_BY_SEQUENCE_HIGH_CONFIDENCE']} associações;
- sucessor funcional de locus fragmentado: {mapping_counts['RESOLVED_SPLIT_FUNCTIONAL_SUCCESSOR']} associações;
- correspondência ambígua: {mapping_counts['AMBIGUOUS_PREVIOUS_STABLE_ID']} associações;
- ainda não encontrado em WBPS19: {mapping_counts['NOT_FOUND_IN_WBPS19']} associações.

Esses números são contagens de associações funcionais, não de evidência biológica.
A evidência de resolução dos cinco IDs está congelada em
`provenance/histone_modifier_id_resolution/historical_id_resolution.tsv`.

## O que não foi misturado ao núcleo

- **HP1:** reader/organizador de heterocromatina, mantido como contexto por decisão do projeto.
- **Histonas canônicas e variantes:** são substratos do sistema, não writers/erasers/readers.
- **Remodeladores e chaperonas:** alteram acesso/organização da cromatina, mas não
  necessariamente escrevem, apagam ou leem uma marca covalente.
- **Genes inferidos apenas por domínio:** somente os dois PHDs inconclusivos
  permanecem em `REVIEW_ADD`; hipóteses contraditas pelos perfis/ortologia foram excluídas.

## Arquivos e ativação

A tabela completa e congelada está em `config/candidates/histone_modifier_gene_set_v2.tsv`.
O arquivo histórico `config/PRJEB32839/candidate_genes.txt` não foi alterado, para
preservar a reprodutibilidade do atlas já construído. A v2 só deve substituir a lista
ativa somente nas análises que declararem explicitamente esta versão. Os cinco IDs
pendentes já foram resolvidos; os dois PHDs de `REVIEW_ADD` continuam exploratórios.

A auditoria dos antigos `REVIEW_ADD`, incluindo HMMER, BLASTP e ortologia, está em
`docs/validacao_review_add_modificadores_histonas.md` e
`provenance/histone_modifier_review_validation/review_add_validation.tsv`.

## Fontes

- Curadoria histórica: lista fornecida e confirmada pela autora do projeto.
- Referência publicada de contexto para a maquinaria de metilação: {CATALOGUE_URL}
- Catálogo publicado de proteínas com bromodomínio: {BROMO_URL}
- Função leitora de H3K4me3 do PHD de TAF3: https://pubmed.ncbi.nlm.nih.gov/18682226/
- Anotação local: `config/reference/gene_annotations_wbps19.tsv` (SM_V10/WBPS19).
"""
    DOC_OUT.write_text(doc, encoding="utf-8")

    provenance = {
        "schema_version": "1.0",
        "gene_set_version": "2.0.0-frozen-p1",
        "status": "FROZEN",
        "frozen_on": "2026-10-09",
        "reference": "Schistosoma mansoni SM_V10 / WormBase ParaSite WBPS19",
        "principle": "Biological curation and identifier mapping status are independent fields.",
        "inputs": {
            SOURCE.relative_to(ROOT).as_posix(): sha256(SOURCE),
            ANNOTATIONS.relative_to(ROOT).as_posix(): sha256(ANNOTATIONS),
            VALIDATION_TABLE.relative_to(ROOT).as_posix(): sha256(VALIDATION_TABLE),
            VALIDATION_MANIFEST.relative_to(ROOT).as_posix(): sha256(VALIDATION_MANIFEST),
            ID_RESOLUTION_TABLE.relative_to(ROOT).as_posix(): sha256(ID_RESOLUTION_TABLE),
            ID_RESOLUTION_MANIFEST.relative_to(ROOT).as_posix(): sha256(ID_RESOLUTION_MANIFEST),
        },
        "outputs": {
            TSV_OUT.relative_to(ROOT).as_posix(): sha256(TSV_OUT),
            DOC_OUT.relative_to(ROOT).as_posix(): sha256(DOC_OUT),
        },
        "counts": {
            "rows": len(rows),
            "unique_input_genes": unique_all,
            "unique_genes_by_scope_tier": unique_by_tier,
            "historical_mapping_associations": dict(mapping_counts),
        },
        "sources": [CATALOGUE_URL, BROMO_URL],
    }
    PROV_OUT.parent.mkdir(parents=True, exist_ok=True)
    PROV_OUT.write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(provenance["counts"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
