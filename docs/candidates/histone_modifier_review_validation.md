# Validation of `REVIEW_ADD` histone-modifier candidates

This is the current English validation record. The original Portuguese record
is retained in [the historical preparation commit](https://github.com/GMiguelAlves/HelixForge-SMansoni/blob/d89067a0197d278dba2ce837c44b7064bc2a76eb/docs/validacao_review_add_modificadores_histonas.md),
not in the current tree. Its original provenance reference is preserved as a
historical assertion. The authoritative evidence table is
`provenance/histone_modifier_review_validation/review_add_validation.tsv`.

## Outcome

The 11 candidates originally classified as `REVIEW_ADD` were evaluated
against the *Schistosoma mansoni* SM_V10/WBPS19 proteome and orthology table,
UniProtKB, Pfam/HMMER profiles, and BLASTP against reviewed human proteins.

| Decision | Genes | Consequence |
|---|---|---|
| Promote to core | `Smp_084140`, `Smp_127180` | Move to `CORE_ADD` |
| Keep exploratory | `Smp_185180`, `Smp_316920` | Remain in `REVIEW_ADD` |
| Exclude from direct panel | `Smp_024380`, `Smp_059920`, `Smp_064250`, `Smp_128520`, `Smp_128530`, `Smp_204920`, `Smp_329500` | Record in `EXCLUDED_REVIEW`; also retain `Smp_064250` as transcriptional context |

## Promotion criterion

A candidate was promoted only when evidence combined a compatible domain,
functional architecture or orthology, and a known function involving histones
or histone marks. A PHD or Tudor domain alone was not sufficient: these
domains also occur in proteins that interact with DNA, RNA, or other proteins.

## Promoted genes

### `Smp_084140` — sirtuin eraser

- 321-aa protein.
- Complete, highly significant SIR2/PF02146 domain.
- UniProt annotation as a class I sirtuin.
- Best human hit: SIRT2, 46.1% identity, 89% coverage,
  `E = 5.46 × 10⁻97`.
- Orthologs in other schistosomes are also annotated as NAD-dependent
  deacetylases.

Decision: `CORE_ADD / ERASER / SIRTUIN_DEACETYLASE`.

### `Smp_127180` — TAF3-PHD reader

- 1,543-aa protein with PHD/PF00628.
- *S. haematobium* ortholog annotated as TFIID subunit 3.
- Best human hit: TAF3.
- TAF3-PHD recognition of H3K4me3 is supported by structural and affinity
  experiments.

Decision: `CORE_ADD / READER / TAF3_PHD_H3K4ME3_READER`.

## Genes retained for review

### `Smp_185180`

The PHD domain is confirmed, but functional orthology is insufficient to
assign recognition of H3K4me3, H3K4me0, or another mark. It remains
exploratory.

### `Smp_316920`

This 161-aa model contains two PHD domains and aligns to PHD modules from
KAT6A/KAT6B and DPF. However, it encodes almost only the tandem module and
lacks reliable whole-protein orthology. It remains exploratory pending gene-
model review or experimental validation of binding specificity.

## Genes excluded from the direct panel

- `Smp_024380`: RMT2 also occurs in GAMT and methyltransferases with
  non-histone substrates. The best human hit was GAMT, not a histone PRMT.
- `Smp_059920`: Only 210 aa and lacks SET/PF00856; insufficient support for
  classification as an EZH/H3K27 methyltransferase.
- `Smp_064250`: Contains PHD, but dominant orthology is TCF20. It may remain
  in transcriptional-regulation context, not as a direct reader.
- `Smp_128520` and `Smp_128530`: No PF02373/JmjC, catalytic HxD motif, or
  curated ortholog compatible with a histone demethylase.
- `Smp_204920` and `Smp_329500`: Neither a Tudor/PF00567 profile nor
  functional orthology supporting histone reading; `Smp_329500` is a 70-aa
  model.

## Methods and parameters

- Reference: SM_V10, WormBase ParaSite WBPS19.
- Proteome: `schistosoma_mansoni.PRJEA36577.WBPS19.protein.fa.gz`.
- Orthology: `schistosoma_mansoni.PRJEA36577.WBPS19.orthologs.tsv.gz`.
- HMMER 3.4 with Pfam profiles PF00567, PF00628, PF00856, PF02146,
  PF02373, PF12796, PF13771, and PF27704.
- BLASTP 2.12.0+ against human UniProtKB/Swiss-Prot, `E ≤ 1e-5`, up to ten
  hits.
- All 11 loci had one protein transcript in the WBPS19 FASTA used.

## External sources

- WormBase ParaSite WBPS19: https://parasite.wormbase.org/ftp.html
- UniProtKB REST: https://rest.uniprot.org/
- InterPro/Pfam: https://www.ebi.ac.uk/interpro/
- TAF3-PHD recognition of H3K4me3:
  https://pubmed.ncbi.nlm.nih.gov/18682226/
