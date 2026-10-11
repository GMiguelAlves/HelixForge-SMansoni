# PRJNA312093 results

The native single-end RNA-seq workflow completed successfully with HelixForge
`v1.1.0`: 161/161 tasks finished, all 20 ENA FASTQs passed integrity checks,
and all 20 samples passed the prespecified technical QC thresholds. Salmon,
tximport, five frozen DESeq2 contrasts, MultiQC and the terminal RNA-seq
manifest were produced against SM_V10/WBPS19.

**Scientific status: `PASS_WITH_LIMITATIONS` for the four non-adult contrasts.**
The two samples labelled adult male correlate at only 0.800 in VST space,
whereas the other nine replicate pairs correlate at 0.943–0.994. Both adult
male samples are closer to samples labelled cercaria than to each other.
Input checksums and sample-to-run mapping passed, so this is not evidence of a
simple file swap within HelixForge; the original biological labels or another
source of discordance still require investigation. No samples were relabelled
or excluded from the original run. The adult female-versus-male contrast
remains `NOT_EVALUABLE`.

The frozen single-end fragment-length prior was 200/80. A prespecified
alternative of 250/80 was run in a separate private sensitivity namespace;
the five contrasts were stable, but that does not resolve sample identity.
The optional sensitivity archive and complete re-entry evidence remain in
private project storage. A separate DESeq2 refit using the 16 non-adult
samples retained the frozen statistical settings. Its eight replicate pairs
correlate at 0.945–0.994 (VST); shared significant genes do not reverse
direction. Significant-gene membership changes, especially for stages 1 and
3, so this is accepted with limitations rather than as an unchanged result.
The compact public package contains:

- `qc/`: standardized per-sample and study-level QC. Historical fields named
  `*_both_mates` count the available R1 reads; there is no R2.
- `expression/`: gene-level counts, TPM and effective-length matrices plus
  import metadata.
- `differential_expression/`: original **provisional** 20-sample contrast
  tables and figures, retained for audit, not the accepted inferential source.
- `reanalysis/nonadult_16/`: the four accepted-with-limitations non-adult
  contrast tables, sample inclusion/exclusion records, replicate correlations,
  PCA coordinates, and comparison against the original fit.
- `manifests/`: terminal and reference provenance. The terminal validation
  was performed against the complete private result tree, not this compact
  publication view.
- `reports/`: sanitized MultiQC, execution review and gate record.

FASTQs, FastQC HTML/ZIP files, complete Salmon directories, RDS files, cache,
workdirs and scheduler logs are excluded. The public manifest is provenance;
re-entry relies on the separately preserved private archive.

See [the execution review](reports/PRJNA312093_execution_report.md) and
[the refit contract](reanalysis/nonadult_16/README.md). The atlas has not yet
been rewired to the new tables; until then it must not ingest the original
20-sample differential-expression tables as accepted evidence.

`PRJNA312093_TECHNICAL_EXECUTION = PASS`

`PRJNA312093_NONADULT = PASS_WITH_LIMITATIONS`

`PRJNA312093_ADULT_CONTRAST = NOT_EVALUABLE`
