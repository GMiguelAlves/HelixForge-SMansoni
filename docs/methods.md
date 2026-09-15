# Methods

## Current RNA-seq methodology

This project executes the released HelixForge v1.0.1 RNA-seq workflow. The
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

## Current boundary

Repository initialization and metadata/reference curation are complete. The
first execution phase begins with RNA-seq. ChIP-seq, RNA-seq/ChIP-seq
integration, cross-study meta-analysis, atlas construction, coexpression, and
candidate prioritization remain planned parts of the repository, but require
separate preregistration before execution.
