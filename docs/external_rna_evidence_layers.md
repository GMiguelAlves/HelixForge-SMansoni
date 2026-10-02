# External RNA evidence and single-end study policy

This document freezes the role of three public RNA studies considered for the
*Schistosoma mansoni* atlas. Only Picard et al. is a new bulk atlas study.
Anderson et al. and Wang et al. remain external evidence layers and do not add
samples, TPM columns, PCA coordinates or DESeq2 contrasts to the bulk atlas.

| Study | Project role | Raw reads processed by HelixForge bulk RNA-seq? |
| --- | --- | --- |
| Picard et al. — `PRJNA312093` | `BULK_ATLAS_STUDY` | Yes, through the isolated single-end Salmon path |
| Anderson et al. — `PRJNA294789` | `ANNOTATION_SUPPORT_ONLY` | No |
| Wang et al. — `PRJNA395457` | `CELL_TYPE_MARKER_REFERENCE` | No |

## PRJNA312093 — Picard et al.

### Atlas role

`BULK_ATLAS_STUDY`. The study contributes the missing developmental design:
male and female worms across five stages. It can support within-stage
female-versus-male contrasts and trajectories of sex-biased expression, while
remaining a separate study rather than being merged with paired-end libraries.

### Processing policy

The raw reads will be processed against SM_V10/WBPS19 through the released
HelixForge single-end path:

`FASTQ → FastQC → Trim Galore → technical-run merge → Salmon → tximport → DESeq2 → reports`

Before execution, the project must freeze and validate:

1. run-to-biological-sample mapping and biological replicate identities;
2. stage and sex labels without collapsing biologically distinct groups;
3. one single-end layout for the complete dataset;
4. fragment-length mean and standard deviation, with their source documented;
5. contrasts within each stage and an estimable DESeq2 design;
6. the HelixForge release/commit and SM_V10/WBPS19 reference checksums.

The study enters the atlas only after its QC, manifests, DE outputs and
provenance pass the project acceptance checks. Slurm certification of the new
HelixForge path is an operational gate and is intentionally separate from this
study-selection decision.

Reference: [Picard et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC5038963/)

## Common external-evidence contract

Every imported record must retain the source accession, publication, source
table or figure, original identifier, mapping method, mapping-source version,
mapping status, source checksum and exact assertion supported by the source.
The canonical join key is an SM_V10/WBPS19 gene identifier.

Historical identifiers may be converted automatically only when the declared
mapping source gives one unique previous-stable-ID mapping. Ambiguous and
unmapped records remain visible as such; they are not assigned to the nearest
or most plausible gene. Published differential-expression labels remain claims
of the source study and are never converted into HelixForge DESeq2 results.

The first curated release should produce:

- `external_annotation_evidence.tsv`: `source_study`, `source_publication`,
  `source_location`, `evidence_type`, `original_id`, `current_gene_id`,
  `mapping_method`, `mapping_source`, `mapping_source_version`,
  `mapping_status`, `stage_or_context`, `assertion`, `source_checksum`,
  `curation_status`, and `notes`;
- `cell_type_marker_evidence.tsv`: `source_study`, `source_publication`,
  `source_location`, `original_id`, `current_gene_id`, `mapping_method`,
  `mapping_source`, `mapping_source_version`, `mapping_status`,
  `cell_population`, `lineage_or_state`, `life_stage`, `marker_class`,
  `reported_statistic`, `assertion`, `source_checksum`, `curation_status`,
  and `notes`.

## PRJNA294789 — Anderson et al.

### Atlas role

`ANNOTATION_SUPPORT_ONLY`. The Roche 454 libraries, one aggregate RNA-seq
library per condition and absence of independent biological replicates make the
study unsuitable for the main quantitative atlas or DESeq2. Its value is the
published transcript reconstruction and historical gene-model evidence.

Accepted material includes gene-model extensions, merge/split proposals,
putative novel genes, TSA accessions and qualitative egg-, female- or
male-enrichment claims. Each assertion must retain its original source and
mapping status. No 454-derived TPM matrix, simulated replicate or cross-study
trajectory is produced.

Reference: [Anderson et al.](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004334)

## PRJNA395457 — Wang et al.

### Atlas role

`CELL_TYPE_MARKER_REFERENCE`. This single-cell study provides marker and cell
context for interpreting bulk candidates; it is not treated as a bulk study
and will not be run through Salmon → tximport → DESeq2 as if cells were
biological replicates.

The curated layer may contain published stem, germinal and differentiated cell
markers, their lineage/state and stage context, ranking/statistic and source
table or figure. The atlas may report that a bulk candidate is a published
marker, but may not infer a cell-abundance change without an independently
reviewed composition model.

Any future raw-data analysis requires a separate single-cell workflow with its
own QC, normalization, embeddings, cell-state annotation, differential model
and provenance.

Reference: [Wang et al.](https://elifesciences.org/articles/35449)

## Acceptance checkpoints

An external-evidence release is accepted only when every row has a publication,
source location, mapping method, mapping-source version, checksum and mapping
status; every populated current identifier belongs to SM_V10/WBPS19; ambiguous
mappings remain unresolved; and no external record changes bulk
study/sample/contrast counts. Picard is excluded from these evidence tables
because it is a primary bulk study with its own terminal manifest.
