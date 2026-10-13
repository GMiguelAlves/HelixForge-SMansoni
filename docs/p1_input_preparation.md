# RNA-seq input preparation: historical P1 batch

## Scope of this preparation record

This record describes preparation of seven accessions that had not yet been
processed when the input package was frozen. It freezes metadata,
configurations, contrast decisions, and provenance; it did **not** execute
these studies or modify the atlas. P1 denotes a preparation batch, not a
separate scientific project. Current per-study status is maintained in the
[study registry](studies_and_preparation.md).

| Accession | Runs | Biological samples | Layout | Mode | Prespecified DE decision |
|---|---:|---:|---|---|---|
| `PRJNA312093` | 20 | 20 | single | `full` | Five female-versus-male contrasts within stage |
| `PRJNA494090` | 24 | 24 | paired | `full` | Selected versus control at P6, P7, and P8 |
| `PRJNA704646` | 24 | 8 | paired | `import` | Blocked: three technical libraries from each RNA pool |
| `PRJNA1150591` | 15 | 15 | paired | `full` | MCLZ versus control in seven-week adults; the four-week cohort was reserved for later linkage with `PRJNA597909` |
| `PRJNA1177852` | 8 | 8 | paired | `import` | Blocked: one library per time point, without parallel controls |
| `PRJNA551487` | 10 | 10 | paired | `full` | FPL-64176 versus DMSO, conditional on the QC gate |
| `PRJNA428470` | 6 | 6 | paired | `full` | 5-azacitidine versus control in adult females |

The inventory totaled 107 runs, 91 biological units, and 352,540,122,551
bytes of raw FASTQ according to ENA snapshots frozen on October 9, 2026.

The meclonazepam paper states that its five four-week MCLZ samples should be
paired with the already processed four-week cohort deposited under
`PRJNA597909`. That cross-deposit contrast was not executed during input
preparation; it requires a later integration analysis that preserves both
studies' provenance.

## Input contracts

- `metadata/rnaseq/<ACCESSION>/runs.tsv`: URL, MD5, size, and run-to-sample
  relationship.
- `metadata/rnaseq/<ACCESSION>/samples.tsv`: biological unit for inference.
- `metadata/rnaseq/<ACCESSION>/metadata.csv`: native HelixForge input.
- `config/<ACCESSION>/de_spec.json`: permitted contrasts.
- `config/<ACCESSION>/de_spec.blocked.json`: explicit inference block.
- `provenance/<ACCESSION>/source_snapshots/ena_read_run.tsv`: raw ENA response.
- `provenance/<ACCESSION>/input_manifest.json`: hashes, selection, and raw-data
  volume.

Derived artifacts and source snapshots are frozen. There is no active P1
generator in this repository; regeneration requires explicit new curation.
`tests/test_p1_preparation.py` checks counts, layouts, download contracts,
contrast replication, and source hashes.

## Single-end input

`PRJNA312093` declares `library_layout=single`,
`fragment_length_mean=200`, and `fragment_length_sd=80` as required by the
then-current HelixForge contract. These values are a provisional protocol
prior, not a measurement reported by the paper. Biological execution must
assess sensitivity to this prior before accepting the study into the atlas.
The STAR path was not selected; the supported path was FastQC → Trim Galore →
merge → Salmon → tximport → DESeq2.

## Frozen gene list

`config/candidates/histone_modifier_gene_set_v2.tsv` preserves historical
curation and adds previously missing classes, separating writers, erasers,
readers, and contextual groups. Five previously unresolved historical IDs
were resolved and audited in
`provenance/histone_modifier_id_resolution/historical_id_resolution.tsv`.
The two PHD proteins in `REVIEW_ADD` remain exploratory, while
`EXCLUDED_REVIEW` records are outside the candidate core.
The [English curation record](histone_modifier_gene_set_v2.md) describes the
gene-set tiers and links to the unchanged frozen evidence.
