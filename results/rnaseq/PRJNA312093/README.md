# PRJNA312093 results

The native single-end RNA-seq workflow completed successfully with HelixForge
`v1.1.0`: 161/161 tasks finished, all 20 ENA FASTQs passed integrity checks,
and all 20 samples passed the prespecified technical QC thresholds. Salmon,
tximport, five frozen DESeq2 contrasts, MultiQC and the terminal RNA-seq
manifest were produced against SM_V10/WBPS19.

**Scientific status: `REVIEW`; study promotion and atlas inclusion: `HOLD`.**
The two samples labelled adult male correlate at only 0.800 in VST space,
whereas the other nine replicate pairs correlate at 0.943–0.994. Both adult
male samples are closer to samples labelled cercaria than to each other.
Input checksums and sample-to-run mapping passed, so this is not evidence of a
simple file swap within HelixForge; the original biological labels or another
source of discordance still require investigation. No samples were relabelled
or excluded. The adult female-versus-male contrast, and any synthesis using
this study, must not be presented as accepted biological evidence yet.

The frozen single-end fragment-length prior was 200/80. A prespecified
alternative of 250/80 was run in a separate private sensitivity namespace;
the five contrasts were stable, but that does not resolve sample identity.
The optional sensitivity archive and complete re-entry evidence remain in
private project storage. The compact public package contains:

- `qc/`: standardized per-sample and study-level QC. Historical fields named
  `*_both_mates` count the available R1 reads; there is no R2.
- `expression/`: gene-level counts, TPM and effective-length matrices plus
  import metadata.
- `differential_expression/`: **provisional** frozen contrast tables and
  figures; these are outputs, not approved biological conclusions.
- `manifests/`: terminal and reference provenance. The terminal validation
  was performed against the complete private result tree, not this compact
  publication view.
- `reports/`: sanitized MultiQC, execution review and gate record.

FASTQs, FastQC HTML/ZIP files, complete Salmon directories, RDS files, cache,
workdirs and scheduler logs are excluded. The public manifest is provenance;
re-entry relies on the separately preserved private archive.

See [the execution review](reports/PRJNA312093_execution_report.md) for the
QC, prior-sensitivity and replicate-concordance evidence.

`PRJNA312093_TECHNICAL_EXECUTION = PASS`

`PRJNA312093_BIOLOGICAL_REVIEW = HOLD`
