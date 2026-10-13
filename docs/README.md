# Documentation

Active project documentation is written in English. The
[study registry](studies_and_preparation.md) records current study status;
the [ChIP-seq core preparation record](chipseq_core_preparation.md) documents
reviewed ChIP-seq selection and execution decisions. The
[RNA-seq input-preparation record](p1_input_preparation.md) describes a
historical preparation batch, not a separate scientific project.

The [histone-modifier gene-set record](histone_modifier_gene_set_v2.md) and
[candidate review](histone_modifier_review_validation.md) are available in
English. Their original Portuguese documents were removed from the current
tree and remain retrievable at commit
[`d89067a`](https://github.com/GMiguelAlves/HelixForge-SMansoni/commit/d89067a0197d278dba2ce837c44b7064bc2a76eb).
Frozen provenance JSON still names those historical paths (and records a
checksum for the gene-set document); those entries were not rewritten to
claim that the English translations are byte-identical. The recorded gene-set
document checksum used CRLF bytes, whereas Git stores the text with LF.

Assay-specific file layouts are described in
[`metadata/rnaseq/`](../metadata/rnaseq/README.md) and
[`metadata/chipseq/`](../metadata/chipseq/README.md).
