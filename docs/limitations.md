# Limitations

- Public studies differ in organism stage, tissue, treatment, library
  preparation, sequencing date, depth, and metadata quality. Cross-study
  heterogeneity must be modeled explicitly before any combined inference.
- Historical protocol differences may be inseparable from biological study
  effects.
- The native HelixForge v1 RNA-seq production path is certified for paired-end
  libraries. Single-end support is a separate HelixForge feature and is not
  assumed here.
- Statistical designs are study-specific. A contrast valid in one project is
  not automatically portable to another.
- PRJNA602528 has one public library per time point. Its frozen execution stops
  after Import; DESeq2 is not estimable under the current categorical design.
- WormBase ParaSite annotation evolves. All reported identifiers and mappings
  must retain the exact WBPS19 annotation identity used here.
- Some public metadata fields are ambiguous or absent. `not_declared` and
  `not_applicable` are preserved instead of being guessed.
- PRJEB32839 is frozen to the 150 runs in the author-curated metadata. Thirty
  additional runs associated with the current ENA project record remain
  explicitly excluded.
- The HelixForge v1 `-resume` behavior was environment-dependent for complete
  workflows on shared HPC storage. Cache reuse must be verified operationally.
- No cross-study meta-analysis, atlas, coexpression, candidate prioritization,
  or epigenomic integration has been validated by this repository yet.

