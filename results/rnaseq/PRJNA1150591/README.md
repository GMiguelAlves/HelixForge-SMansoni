# PRJNA1150591 results

HelixForge `v1.1.0` completed the paired-end RNA-seq path for all 15 biological
libraries: QC, Salmon, Import, DESeq2, MultiQC and terminal manifest. The
reference was SM_V10/WBPS19. All input FASTQs passed integrity checks, all
samples passed the predefined technical QC thresholds, and the workflow
completed 167/167 tasks.

The frozen model is `~ condition` and includes all three recorded groups:
five adult 7-week controls, five adult 7-week meclonazepam-treated samples,
and five juvenile 4-week treated samples. The juvenile samples are **retained
in the model because they belong to this study**; the sole published contrast
is treated versus control **within 7-week adults**. This does not imply a
juvenile-versus-adult treatment effect, and no such contrast was tested.

The adult contrast tested 9,458 genes; 936 met the frozen significance rule.
These are statistical outputs pending biological interpretation, not a claim
that all 936 are validated treatment targets. The untreated 4-week comparison
belongs to a separate future cross-study analysis and is not inferred here.

This compact public package contains standardized QC, gene-level counts/TPM/
length matrices, differential-expression tables and figures, source manifests,
sanitized MultiQC and an [execution review](reports/PRJNA1150591_execution_report.md).
The terminal manifest validation refers to the complete original result tree,
not this reduced publication view. FASTQs, individual FastQC reports, complete
Salmon directories, RDS objects, workdirs, scheduler logs and the private
re-entry audit archive remain outside Git.

`TECHNICAL_EXECUTION = PASS`

`SCIENTIFIC_INTERPRETATION = PENDING_REVIEW`
