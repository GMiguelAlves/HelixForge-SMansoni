# Histone-modifier gene set v2 (SM_V10/WBPS19)

This is the current English curation record. The original Portuguese record
is retained in [the historical preparation commit](https://github.com/GMiguelAlves/HelixForge-SMansoni/blob/d89067a0197d278dba2ce837c44b7064bc2a76eb/docs/conjunto_modificadores_histonas_v2.md),
not in the current tree. Its frozen provenance entry remains a historical
assertion, not a checksum for this translation. The authoritative data table is
`config/candidates/histone_modifier_gene_set_v2.tsv`.

## Curation decision

The historical list is a biological curation of this project. Members of
`Writers`, `Readers`, and `Erasers` are treated as histone modifiers confirmed
by the project author's curation. HP1 is retained separately as a
heterochromatin reader/organizer and is excluded from the enzymatic-modifier
count.

Identifier status is distinct from biological classification. The five
historical IDs that were previously ambiguous or absent were resolved through
sequence comparison and, for fragmented loci, preservation of the expected
functional segment.

## Gene-set tiers

| Tier | Criterion | Unique genes | Atlas role |
|---|---|---:|---|
| `CORE_CURATED` | Project-curated historical list, excluding HP1 from the enzymatic count | 60 | Core after resolution of pending IDs |
| `CORE_ADD` | Explicit annotation, published catalog, or domain-plus-orthology validation | 38 | Main expansion |
| `REVIEW_ADD` | Confirmed domain but unresolved histone function or gene model | 2 | Exploratory panel outside the core |
| `CONTEXT` | HP1, cofactors, remodelers, and chaperones | 24 | Mechanistic context outside the direct-modifier count |
| `EXCLUDED_REVIEW` | Initial hypothesis not supported as a direct modifier/reader | 7 | Audited but not analyzed as a candidate |

There are **126 unique input genes** across all tiers. The table records
functional associations, so a gene may occupy more than one row—for example,
an acetyltransferase that also contains a bromodomain.

## Notable additions and decisions

- **Acetylation writers:** HATs explicitly annotated in SM_V10 but absent
  from the historical list.
- **Acetylation erasers:** Classical HDACs and sirtuins. `Smp_084140` was
  promoted after confirmation of a SIR2 domain and SIRT2 orthology.
- **Acetylation readers:** 22 proteins from the published bromodomain catalog.
- **H3K4me3 reader:** `Smp_127180` was promoted as a PHD-containing TAF3
  ortholog; TAF3-PHD recognition of H3K4me3 is experimentally established.
- **Exploratory PHD proteins:** `Smp_185180` and `Smp_316920` have PHD domains
  but insufficient support for histone-mark specificity.
- **Rejected direct-panel hypotheses:** `Smp_024380`, `Smp_059920`,
  `Smp_064250`, `Smp_128520`, `Smp_128530`, `Smp_204920`, and `Smp_329500`
  remain in `EXCLUDED_REVIEW` for audit.
- **Chromatin context:** Polycomb/Trithorax, CoREST/SAP18/RBBP4,
  ATP-dependent remodelers, and histone chaperones have separate groups.

## Functional-association counts

| Tier | Class | Gene–class associations |
|---|---|---:|
| `CORE_CURATED` | `ERASER` | 13 |
| `CORE_CURATED` | `READER` | 24 |
| `CORE_CURATED` | `WRITER` | 26 |
| `CORE_ADD` | `ERASER` | 9 |
| `CORE_ADD` | `READER` | 23 |
| `CORE_ADD` | `WRITER` | 10 |
| `REVIEW_ADD` | `READER` | 2 |
| `CONTEXT` | `CHAPERONE` | 4 |
| `CONTEXT` | `COFACTOR` | 10 |
| `CONTEXT` | `HP1` | 2 |
| `CONTEXT` | `REMODELER` | 8 |
| `EXCLUDED_REVIEW` | `ERASER` | 2 |
| `EXCLUDED_REVIEW` | `READER` | 3 |
| `EXCLUDED_REVIEW` | `WRITER` | 2 |

## Historical-ID resolution

- Exact SM_V10 match: 48 associations.
- Unique previous stable-ID match: 12 associations.
- Exact sequence match: 2 associations.
- High-confidence sequence match: 1 association.
- Functional successor of a fragmented locus: 2 associations.
- Ambiguous match: 0 associations.
- Not found in WBPS19: 0 associations.

These are functional-association counts, not counts of independent biological
evidence. Resolution evidence for the five IDs is frozen in
`provenance/histone_modifier_id_resolution/historical_id_resolution.tsv`.

## Exclusions from the core

- **HP1:** Heterochromatin reader/organizer retained as context by project
  decision.
- **Canonical and variant histones:** Substrates of the system, not
  writers, erasers, or readers.
- **Remodelers and chaperones:** Alter chromatin access or organization but
  do not necessarily write, erase, or read a covalent mark.
- **Genes inferred from domains alone:** Only the two unresolved PHD proteins
  remain in `REVIEW_ADD`; hypotheses contradicted by profiles or orthology
  were excluded.

## Activation and supporting evidence

The frozen table is `config/candidates/histone_modifier_gene_set_v2.tsv`.
The historical `config/rnaseq/PRJEB32839/candidate_genes.txt` was not changed, so
the previously built atlas remains reproducible. Analyses must explicitly
declare v2 before replacing their active candidate list. The five historical
IDs are resolved; the two `REVIEW_ADD` PHD proteins remain exploratory.

The [English review of `REVIEW_ADD` genes](histone_modifier_review_validation.md)
summarizes HMMER, BLASTP, and orthology evidence. Its frozen evidence table is
`provenance/histone_modifier_review_validation/review_add_validation.tsv`.

## Sources

- Historical curation: list supplied and confirmed by the project author.
- Published context for methylation machinery:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6125318/
- Published bromodomain-protein catalog:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10726355/
- H3K4me3 recognition by TAF3-PHD:
  https://pubmed.ncbi.nlm.nih.gov/18682226/
- Local annotation: `config/reference/gene_annotations_wbps19.tsv`
  (SM_V10/WBPS19).
