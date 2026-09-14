# Methods

## HelixForge methodology

This project executes the released HelixForge v1.0.0 RNA-seq workflow. The
supported production route is run-level metadata validation, reference-bundle
validation, FastQC, Trim Galore, post-trim QC, technical-run consolidation,
Salmon quantification, tximport under the `production_v1` policy, DESeq2 when
the design is estimable, and terminal manifest generation.

HelixForge owns workflow implementation, process resources, scientific API
contracts, terminal schemas, and generic provenance. The exact software release
and commit are pinned in this repository and must be verified before every
run. This project does not copy or modify the HelixForge core.

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

## Current boundary

Repository initialization and metadata/reference curation are complete.
Scientific processing has not started. Cross-study meta-analysis, atlas
construction, coexpression, candidate prioritization, and epigenomic
integration require separate preregistration and are out of scope here.

