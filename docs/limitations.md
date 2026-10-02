# Limitations

- Public studies differ in organism stage, tissue, treatment, library
  preparation, sequencing date, depth, and metadata quality. Cross-study
  heterogeneity must be modeled explicitly before any combined inference.
- Historical protocol differences may be inseparable from biological study
  effects.
- The current accepted studies use the paired-end HelixForge path. Native
  single-end Salmon support is implemented upstream but awaits isolated Slurm
  certification and a released project pin before `PRJNA312093` is executed.
  No existing paired-end study is reclassified or mixed with that dataset.
- `PRJNA294789` is reserved for historical annotation evidence and
  `PRJNA395457` for published cell-type markers. Neither is eligible for the
  bulk Salmon/tximport/DESeq2 route under the current atlas design.
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
- Nextflow cache reuse requires the official launcher, a stable launch/cache
  location, the same work directory, and preserved task outputs. Controlled
  validation with Nextflow 25.10.7 and Java 21 recovered every eligible task
  and confirmed selective invalidation. The earlier failures attributed to
  shared NFS were caused by unsupported direct JAR invocation and are resolved.
- In PRJNA597909, the scientific workflow completed, but terminal-manifest
  generation initially selected a host Python without the required
  `jsonschema` package. Runtime ordering was corrected and only the terminal
  boundary was recovered under Slurm; QC, Salmon, Import, and DESeq2 were not
  recomputed.
- A historical PRJNA597909 `-resume` attempt failed to reuse eligible entries
  and was stopped before expensive stages repeated. Later investigation showed
  that the run bypassed launcher-supplied JVM serialization settings by calling
  the Nextflow JAR directly. This incident remains part of the study audit
  trail, but it is not a current shared-NFS limitation.
- The RNA-seq atlas is descriptive. Its global PCA and correlation can reflect
  both biology and study/protocol effects and must not be interpreted as formal
  batch-effect estimates.
- Batch metadata availability and interpretation are study-specific. Formal
  cross-study batch assessment is deferred and no corrected matrix is
  presented.
- No ChIP-seq application, cross-study inferential meta-analysis, coexpression,
  candidate prioritization, or RNA-seq/ChIP-seq integration has been validated
  by this repository yet.
