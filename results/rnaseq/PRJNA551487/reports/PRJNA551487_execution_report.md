# PRJNA551487 — execution and QC review

## Completed path

The paired-end HelixForge `v1.1.0` run completed 117/117 tasks for ten
biological libraries (five DMSO, five FPL). All 20 FASTQs passed size,
checksum, gzip and pairing validation. Salmon 1.10.3, tximport, MultiQC and
the frozen DESeq2 contrast executed with the SM_V10/WBPS19 reference. The
count, TPM and effective-length matrices contain 9,914 genes × ten
libraries; their checksums match the native
[Import manifest](../expression/import_manifest.json).

## Independent QC gate

| Measure | Observed |
| --- | ---: |
| Libraries classified `PASS` | 0/10 |
| Libraries classified `REVIEW` | 10/10 |
| Post-trim retention | 99.75–99.81% |
| Salmon mapping | 12.41–15.02% |
| Salmon mapping mean | 13.96% |
| Mapped/processed fragments | 23,514,623/168,740,503 |
| MultiQC | `PASS` |

The predefined review triggers are post-trim retention below 80% and Salmon
mapping below 60%. Only the mapping trigger fired. The mapped fraction was
12.41–13.81% in DMSO and 14.70–15.02% in FPL. The
[exploratory PCA](replicate_pca.png) and
[replicate correlations](replicate_correlations.tsv) describe the mapped
counts, but cannot rule out bias in the much larger unmapped fraction.

A preliminary BLAST-based diagnostic sampled reads only from the beginnings
of the FASTQs, so its hit fractions were **not** treated as whole-library
estimates. Before inspecting the subsequent results, a private protocol
froze a deterministic, uniform reservoir of 3,000 paired reads from each of
the ten complete FASTQ pairs. Each original source passed its frozen ENA size and
MD5 checks, mates remained paired, and every positional decile contributed
to each sample (at least 264 and at most 337 of the 3,000 pairs per decile).
Two initial transfers failed MD5 despite reaching the expected byte size;
only those transfers were repeated from byte zero and accepted after the
checksum matched. The other eight samples were not repeated.

In this 30,000-pair diagnostic, 86.60–92.87% of pairs had at least one mate
locally aligned by Bowtie2 (`--very-sensitive-local`, without a mapping-quality
filter) to the predeclared 50-kb repetitive locus
`SM_V10_2:38300000-38350000`. The mean was 91.91% in controls and 88.27% in treated
libraries; these are **screening fractions, not an inferential contrast**.
The public [per-library summary](fullfile_locus_qc_summary.tsv) reports the
denominators and hits. This target was selected after the preliminary
diagnostic, and permissive local alignment to one repetitive locus does not
establish unique origin, rRNA identity, or a specific laboratory cause.
The complete sampled-read, alignment, checksum, tool-version, and execution
audit is preserved privately, not distributed with the compact result package.

## Scientific decision

The atlas decision was recorded without inspecting or selecting favorable
DESeq2 hits. The DESeq2 outputs are **not promoted as biological results** and
remain in private audit storage. No contrast table is included here, no
sample has been excluded, and no threshold was changed. In the
[primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6796685/), DMSO
controls came from mesenteries and FPL-treated worms from livers. Treatment
and collection site therefore cannot be separated by this design. The
study is `EXCLUDE_FROM_ATLAS` for both expression and differential-effect
displays. This is a limitation of comparability and interpretation, **not**
a failure of the workflow or proof that the original publication is wrong.

The [terminal manifest](../manifests/rnaseq_run_manifest.json) and its
[validation](../manifests/terminal_manifest_validation.json) refer to the
complete private execution; the compact public package intentionally omits
DE artifacts after the scientific review withheld their promotion.
