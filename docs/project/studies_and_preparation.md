# Study registry and preparation

This is one HelixForge-SMansoni scientific project. The historical `P1`
label identifies an execution/preparation batch, not a separate project or
biological analysis. Study-specific results enter shared analyses only after
individual review; distinct studies are not automatically pooled into one
DESeq2 model.

## Runtime provenance

The earlier completed studies retain their individual `v1.0.1`/`v1.0.2`
execution provenance. The historical project default was HelixForge
[`v1.0.2`](https://github.com/GMiguelAlves/HelixForge/releases/tag/v1.0.2),
commit `5d4b3e696319db5cd7633472504964f1dc7c0434`.
The later execution template pinned HelixForge
[`v1.1.0`](https://github.com/GMiguelAlves/HelixForge/tree/v1.1.0),
commit `062001ceb46d8a25c2415823b97d0e310658f8c0`, including the native
single-end path. `PRJNA704646` additionally used the reviewed Import fix at
commit `9fb09ff`; its per-study record is authoritative. These are runtime
pins for particular executions, not distinct scientific projects. Nextflow
`25.10.7` and Java `21` are the certified runtime baseline.

## RNA-seq studies

The table is a single study inventory, not separate scientific cohorts.
"Prepared" means that inputs were curated; it does not mean the workflow ran
or that results were accepted for the atlas. Full evidence and restrictions
remain in each published study package.

| Study | Runs | Biological samples | Published/review status |
|---|---:|---:|---|
| `E-ERAD-478` | 120 | 60 | Full RNA-seq; `PASS` |
| `E-MTAB-451` | 12 | 11 | Full RNA-seq; `PASS` |
| `PRJEB14695` | 138 | 23 | Full RNA-seq; `PASS_WITH_LIMITATIONS` |
| `PRJEB3190` | 20 | 20 | Full RNA-seq; `PASS` |
| `PRJEB32839` | 150 | 75 | Full RNA-seq; `PASS_WITH_LIMITATIONS` |
| `PRJNA1150591` | 15 | 15 | Full RNA-seq; adult contrast admitted to atlas with a design warning |
| `PRJNA1177852` | 8 | 8 | Import completed; admitted as a descriptive time series with warning; no DE |
| `PRJNA312093` | 20 | 20 | Technical execution passed; three restricted-refit contrasts admitted with warning, adult/cercaria not evaluable |
| `PRJNA428470` | 6 | 6 | Full RNA-seq completed; coding-gene contrast admitted with mapping warning |
| `PRJNA494090` | 24 | 24 | Technical execution passed; excluded from atlas-level DE interpretation |
| `PRJNA551487` | 10 | 10 | Technical execution passed; excluded from atlas after whole-file mapping investigation and design review |
| `PRJNA597909` | 20 | 20 | Full RNA-seq; `PASS_WITH_LIMITATIONS` |
| `PRJNA602528` | 10 | 10 | Import completed; `PASS_WITH_LIMITATIONS`; DE not applicable by design |
| `PRJNA704646` | 24 | 8 | Import completed; descriptive atlas inclusion with QC warning for three pools |

The frozen input package for seven of these accessions retains its historical
batch label in [`metadata/rnaseq/p1_pending_manifest.tsv`](../../metadata/rnaseq/p1_pending_manifest.tsv)
and [the input-preparation record](../rnaseq/p1_input_preparation.md). It records runs,
biological units, configurations, DE decisions and ENA source snapshots;
preparation alone does not determine analysis or atlas inclusion.

`PRJNA312093` completed the native single-end path under HelixForge `v1.1.0`:
161/161 tasks and 20/20 technical QC checks passed, and the prespecified
fragment-prior sensitivity check passed. Later independent-marker checks found
adult- and cercaria-labelled profile inconsistencies. Only three
schistosomulum-stage sex contrasts from a 12-sample restricted refit are
approved for atlas inclusion, with explicit warnings; its original 20-sample
DE tables are **not** that evidence. Adult and cercaria contrasts remain not
evaluable. See the [study admission decision](../rnaseq/atlas_study_decisions.md)
for the investigation and limits. The verified
[restricted-refit tables](../../results/rnaseq/PRJNA312093/restricted_refit/)
are the only approved DE source; the original 20-sample compact
[result package](../../results/rnaseq/PRJNA312093/) remains an execution and
audit record. Other accessions retain their own review gates.

`PRJNA494090` completed the paired-end workflow, but 10 of 24 libraries fell
below the prespecified 60% Salmon mapping review threshold. Systematic
selected-versus-control differences in quantifiable reads and a variable
external sex-expression signature make atlas-level PZQ-effect interpretation
unsafe. It is excluded from the atlas, not erased from the project; no sample
was automatically excluded. See the [published QC review](../../results/rnaseq/PRJNA494090/README.md)
and [admission decision](../rnaseq/atlas_study_decisions.md).

`PRJNA704646` completed the paired-end Import path for 24 technical runs
assigned to eight biological pools. It produced count, TPM and effective-length
matrices for 9,914 genes. Three pools require QC review; there is one
biological unit per resistance × stage × sex group, so no DESeq2 inference is
claimed. It is admitted to the atlas for descriptive expression only, with a
visible QC warning. See the [published descriptive result](../../results/rnaseq/PRJNA704646/README.md).

`PRJNA1150591` completed the full RNA-seq workflow. Only its prespecified
seven-week adult meclonazepam-versus-control contrast is admitted to the atlas,
with a warning that the five four-week treated juveniles remained in the
original `~ condition` fit. This is not an age contrast. See the
[result package](../../results/rnaseq/PRJNA1150591/README.md) and the
[study admission decision](../rnaseq/atlas_study_decisions.md).

`PRJNA1177852` completed the native Import path for eight post-meclonazepam
time points. All eight libraries passed technical QC, but each time point has
one biological library and no parallel control. It is admitted to the atlas
only as a labelled descriptive trajectory, with no DE or treatment-effect
claim. See the [result package](../../results/rnaseq/PRJNA1177852/README.md)
and [admission policy](../rnaseq/atlas_study_decisions.md).

`PRJNA551487` completed technically, but all ten libraries have Salmon mapping
below the predefined review threshold (12.41%–15.02%). A frozen whole-file
random-pair screen and the published treatment/collection-site design were
reviewed. Its DE outputs remain private and unpromoted, and neither its
expression nor its DE effects enter the atlas. The screen does not establish
rRNA identity or a library-preparation cause. See the
[published QC review](../../results/rnaseq/PRJNA551487/README.md) and
[study decision](../rnaseq/atlas_study_decisions.md).

`PRJNA428470` completed the full paired-end workflow for three adult-female
5-azacitidine libraries and three controls. Four of six libraries remain QC
`REVIEW` for mapping below 60%. A two-library mapping diagnostic and six
leave-one-out sensitivity fits support admitting the frozen contrast within
the **coding-gene mRNA-index scope**, with a visible mapping warning and no
claims about unquantified reads or lncRNAs. See the compact
[result package](../../results/rnaseq/PRJNA428470/README.md) and
[admission policy](../rnaseq/atlas_study_decisions.md).

Every published study result uses the same compact QC interface at
`results/rnaseq/<study>/qc/qc_sample_summary.tsv` and `qc_summary.json`; detailed
MultiQC reports remain available alongside it where published. The fields,
classification rules, and generation procedure are defined in
[QC results contract](../rnaseq/qc_results_contract.md).

`PRJNA602528` is the first completed real *S. mansoni* application and
operational calibration run. It is not a new benchmark of HelixForge. Its
public design has one library per time point, so the frozen analysis stops after
tximport and does not perform DESeq2 inference. The compact accepted products
and execution report are available under [`results/rnaseq/PRJNA602528/`](../../results/rnaseq/PRJNA602528/).

`PRJNA597909` is the first complete native RNA-seq application in this
project. All 20 paired-end samples passed QC, Salmon quantification and
tximport completed, and all four frozen DESeq2 contrasts were estimated. The
accepted tables, figures, portable terminal manifest, and report are available
under [`results/rnaseq/PRJNA597909/`](../../results/rnaseq/PRJNA597909/). Its
`PASS_WITH_LIMITATIONS` classification records terminal-runtime recovery and
the then-unresolved `-resume` incident; the scientific stages themselves
passed. Later investigation traced that incident to unsupported direct JAR
invocation rather than shared NFS, and the official launcher subsequently
passed complete cache-reuse and selective-invalidation validation.

`PRJEB14695` completed the native RNA-seq path for 138 paired-end technical
runs representing 23 biological samples. Run-to-sample mapping and aggregation
were enforced as scientific gates before tximport and DESeq2. The first attempt
exposed an overly strict validation of the run-level `technical_unit` field;
the reviewed HelixForge correction was then validated by a controlled resume
that reused all 815 eligible upstream tasks. The accepted matrices, six frozen
contrasts, figures, HTML reports, terminal manifest, and validation summary are
available under [`results/rnaseq/PRJEB14695/`](../../results/rnaseq/PRJEB14695/).

Earlier studies retain their recorded `v1.0.1` execution provenance.
`PRJEB14695` additionally pins the reviewed post-release contract fix at commit
`42864266892d1165477bb3b33c919e1fabb28ad1`; later studies used `v1.0.2`
without rewriting prior runs.

`PRJEB32839` completed the `v1.0.2` native path for 150 technical runs and 75
biological samples. The exact frozen inventory, 12 DESeq2 contrasts, compact
QC/DE summaries, candidate-gene report record, terminal manifest, storage
calibration, and reviewed execution report are available under
[`results/rnaseq/PRJEB32839/`](../../results/rnaseq/PRJEB32839/). The newer pin applies to this
study and does not rewrite earlier execution provenance. Verified persistence
and the private audit archive preceded cleanup; 1.208 TB of study-specific
scratch data were removed without touching shared resources or other studies.

`E-MTAB-451` completed the native `v1.0.2` paired-end path for 12 runs
representing 11 biological samples. The two technical runs assigned to
`schistosomulum_3h_r1` were consolidated before quantification and inference.
All 136 processes completed, the three frozen `~ condition` contrasts were
estimated, and Salmon mapping rates ranged from 62.26% to 78.87% (median
75.18%). The accepted matrices, differential-expression results, QC summary,
terminal manifests, and sanitized execution reports are available under
[`results/rnaseq/E-MTAB-451/`](../../results/rnaseq/E-MTAB-451/). The complete candidate-gene
report passed but is retained outside Git under the current repository-size
policy. Persistence, checksums, private audit capture, and sensitive-path
scanning passed before the study-specific scratch data and FASTQs were removed.

`PRJEB3190` completed the native `v1.0.2` paired-end path for 20 biological
samples spanning seven schistosomulum time points. All 227 processes completed,
Salmon mapping rates ranged from 63.14% to 88.95% (median 82.95%), tximport
produced matrices for 9,914 genes, and all nine frozen `~ condition` contrasts
were estimated. ENA alias-derived sequencing labels remain exploratory because
they are not documented as technical batches and were not included in the
inferential model. The accepted matrices, differential-expression results, QC
summary, terminal manifests, and sanitized execution reports are available
under [`results/rnaseq/PRJEB3190/`](../../results/rnaseq/PRJEB3190/). The complete candidate-gene
report passed but remains outside Git under the repository-size policy.
Persistence, checksums, private audit capture, and sensitive-path scanning
passed before the study-specific workdir, FASTQs, and other reproducible
scratch intermediates were removed.

`E-ERAD-478` completed the native `v1.0.2` paired-end path for 120 runs
representing 60 biological samples. Two technical runs were aggregated per
sample before quantification and inference. All 950 processes completed,
Salmon mapping rates ranged from 75.03% to 91.84% (median 85.74%), tximport
produced matrices for 9,914 genes, and all 32 frozen contrasts were estimated
with the inferential design `~ batch + condition`. All 60 samples passed the
frozen QC review thresholds. The accepted matrices, differential-expression
results, QC summaries, terminal manifests, and sanitized execution reports are
available under [`results/rnaseq/E-ERAD-478/`](../../results/rnaseq/E-ERAD-478/). The complete
candidate-gene report passed but remains outside Git. Persistence, checksums,
private audit capture, and sensitive-path scanning passed before 621.87 GB of
study-specific scratch data were removed.

## ChIP-seq core preparation

The ChIP-seq study inventory is curated in this same project with
assay-specific control relationships. Five Q1 studies have controls suitable
for peak calling
(`PRJNA1173904`, `PRJNA236156`, `PRJNA312093`, `PRJNA602708`, and
`PRJNA1090249`). `PRJNA219440` is retained for QC/alignment only because
its public historical subexperiments do not provide a uniform matched-control
design. `PRJNA384372` is P2/Q2, and `PRJNA816041` is a technical benchmark
rather than atlas evidence.

The metadata, download checksums, explicit IP/control mapping, narrow/broad
run sets, duplicate-deposit exclusions, and publication links are indexed in
[`metadata/chipseq/`](../../metadata/chipseq/). The frozen execution and reference
decisions are under [`config/chipseq/`](../../config/chipseq/), with the scientific
preparation record in
[ChIP-seq core preparation record](../chipseq/chipseq_core_preparation.md).
No ChIP-seq result has yet been promoted to the atlas.
