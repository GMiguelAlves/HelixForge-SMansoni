# PRJNA428470 results — coding-gene scope

HelixForge `v1.1.0` completed the paired-end RNA-seq path for six biological
libraries (three adult females treated with 5-azacitidine and three controls)
against SM_V10/WBPS19. All 12 FASTQs used in the run passed source checks;
77/77 Nextflow tasks, Salmon, Import, DESeq2, MultiQC and terminal-manifest
validation completed.

The scientific classification is **`PASS_WITH_LIMITATIONS`**, not an
unqualified QC pass. Salmon mapping to the declared **mRNA** index was
34.61%–60.83%; four of six libraries remain `REVIEW` under the frozen 60%
threshold. The treated and control means were 47.08% and 50.11%. A
two-library diagnostic found a larger ribosomal/repetitive and highly
multimapping signal in the lower-mapping library. It does not characterize
all six unmapped fractions or establish contamination by another species.

The prespecified `~ condition` comparison tested 9,307 genes; 1,454 met the
frozen significance criteria. Six supplementary leave-one-out fits retained
the direction of every originally significant gene present in each refit;
1,286–1,363 of the 1,454 remained significant. These fits assess stability
and do **not** replace the full six-library comparison. The gene-level
results may enter the atlas only as coding-gene evidence with a visible
mapping warning. They are not validated drug targets, and the package cannot
support inferences about lncRNAs or other transcripts missing from the index.

The compact package contains [per-sample QC](qc/qc_sample_summary.tsv),
[leave-one-out sensitivity](qc/leave_one_out_summary.tsv),
[gene-level matrices](expression/), the
[complete prespecified contrast](differential_expression/contrasts/DEG_condition__female_adult_5azac_vs_female_adult_control.tsv),
figures, sanitized MultiQC, and portable provenance. See the
[execution and interpretation review](reports/PRJNA428470_execution_report.md).
The terminal-manifest validation applies to the complete original result
tree, not to this reduced publication package. FASTQs, full Salmon
directories, RDS objects, workdirs, scheduler logs and the private audit
archives are excluded from Git.

```ini
TECHNICAL_EXECUTION = PASS
SAMPLE_QC = REVIEW_4_OF_6
DE_CODING_GENE_SCOPE = PASS_WITH_LIMITATIONS
ATLAS_INCLUSION = INCLUDE_WITH_MAPPING_WARNING
```
