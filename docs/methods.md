# Methods

## Current RNA-seq methodology

This project executes the released HelixForge v1.0.2 RNA-seq workflow. The
supported production route is run-level metadata validation, reference-bundle
validation, FastQC, Trim Galore, post-trim QC, technical-run consolidation,
validation and reuse of a prebuilt Salmon index, Salmon quantification,
tximport under the `production_v1` policy, DESeq2 when
the design is estimable, and terminal manifest generation.

HelixForge owns workflow implementation, process resources, scientific API
contracts, terminal schemas, and generic provenance. The exact software release
and commit are pinned in this repository and must be verified before every
run. This project does not copy or modify the HelixForge core.

The preparation preflight originally certified HelixForge `v1.0.0`. Before
the first scientific execution, it was found that the quantification graph
could only reuse an index indirectly through task cache. HelixForge `v1.0.1`
added an explicit manifest-validated prebuilt-index path. The project pin was
updated before scientific execution; no scientific defaults or results were
changed.

## Schistosoma-specific decisions

- Organism: *Schistosoma mansoni* (NCBI Taxonomy 6183).
- Assembly: SM_V10, BioProject PRJEA36577.
- Reference and annotation source: WormBase ParaSite WBPS19.
- Quantification reference: WBPS19 mRNA transcript sequences.
- Annotation inputs: WBPS19 canonical GTF and complete GFF3.
- Import policy: full-length libraries with `countsFromAbundance=lengthScaledTPM`.
- Transcript identifier version/bar stripping remains disabled under the
  HelixForge `production_v1` contract.
- Batch correction is not applied to inferential count matrices. A supported
  batch variable must enter an estimable DESeq2 model formula.

Study-specific conditions, contrasts, exclusions, and design limitations are
kept under `config/<study>/` and `metadata/<study>/`. They are not promoted to
generic HelixForge defaults.

## ChIP-seq and integration methodology

This repository also owns the future *S. mansoni* ChIP-seq applications and
the integration of their accepted evidence with RNA-seq results. Those
analyses will use the released HelixForge ChIP-seq and Integrative workflows,
but will receive independent study metadata, reference compatibility checks,
scientific contracts, provenance, and acceptance reports here.

Integration will consume accepted terminal manifests rather than untracked
paths into upstream work directories. No ChIP-seq or integrative result is
implied by completion of an RNA-seq study.

## RNA-seq atlas methodology

The atlas consumes the accepted sample metadata, TPM matrices, QC summaries,
and original DESeq2 results of all seven completed RNA-seq studies. Compatibility
is gated on the common `Schistosoma_mansoni_SM_V10_WBPS19` identity, identical
gene order, canonical gene identifiers, and biological-sample columns.

Cross-study expression displays use `log2(TPM + 1)`. PCA uses up to the 2,000
most variable genes followed by gene-wise standardization and singular-value
decomposition. Pearson sample correlation uses the complete transformed shared
gene universe. These are descriptive views; no pooled differential-expression
model is fitted. Within-study DESeq2 log2 fold changes and adjusted p-values are
copied without re-estimation and retain numerator/denominator orientation.
Significance markers use each frozen study rule: adjusted p-value below 0.05
and absolute log2 fold change of at least 1.

Gene labels are derived from the frozen WBPS19 GFF3 and retain the canonical
`Smp_*` identifier as the primary key. Functional descriptions, biotype,
source accession, previous stable ID, and coordinates are carried as
annotation metadata. Historical candidate IDs are resolved only when WBPS19
provides one unambiguous current-gene mapping; ambiguous mappings are reported
without selecting an alternative.

The report includes an inventory, QC overview, global and within-study PCA,
sample correlation, selected-gene expression and effect heatmaps, a searchable
gene explorer, exportable SVGs, tables, and checksum provenance. Batch metadata
are retained, but formal batch-effect assessment and corrected visualization
matrices are deferred.

## Current boundary

Repository initialization and metadata/reference curation are complete.
`PRJNA602528` completed QC, Salmon quantification, tximport, terminal-manifest
validation, persistence, and cleanup under its `IMPORT_ONLY` contract. Its
temporal figures are descriptive and no DESeq2 inference was performed.
The seven registered RNA-seq studies have completed their accepted project
contracts. The descriptive RNA-seq atlas is now implemented from their
versioned outputs. ChIP-seq, RNA-seq/ChIP-seq integration, cross-study
inferential meta-analysis, coexpression, and candidate prioritization remain
separate reviewed analyses.
