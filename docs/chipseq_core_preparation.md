# *Schistosoma mansoni* ChIP-seq core preparation

## Scope of the preparation record

This record documents public-data preparation before scientific execution. It
did not change HelixForge, run analyses on the server, or promote results to
the atlas. The core uses the frozen `SM_V10/WBPS19` reference and the native
ChIP-seq contract of HelixForge `v1.1.0`
(`062001ceb46d8a25c2415823b97d0e310658f8c0`).

The ENA inventory was rechecked. On October 12, 2026, WBPS19 remained the
latest version exposed by the official WormBase ParaSite directory. The
existing reference, checksums, and validation were reused; the ChIP-seq
package did not introduce another genome copy.

## Prepared studies

| Study | Scope | Raw FASTQ | Preparation decision |
|---|---:|---:|---|
| `PRJNA1173904` | 68 paired-end runs; eight gonadal marks | 130.77 GB | Q1/P1; peak calling with group-specific inputs |
| `PRJNA236156` | 60 single-end runs; four marks in cercariae/adults | 119.01 GB | Q1/P1; peak calling with input UB0 matched to each biological unit |
| `PRJNA312093` | 24 single-end runs; H3K4me3/H3K27me3 | 34.70 GB | Q1/P1; sex and life cycle, two replicates |
| `PRJNA602708` | 24 paired-end runs; H3K4me3/H3K27me3 | 13.73 GB | Q1/P1; plasma/treatment, two replicates |
| `PRJNA1090249` | Eight *S. mansoni* single-end runs; HP1 | 6.19 GB | Specialized Q1/P1; four *Drosophila* runs excluded |
| `PRJNA219440` | 29 mixed-layout runs; historical collection | 49.43 GB | Q2/P1; QC/alignment only |
| `PRJNA384372` | Three single-end H3K9ac runs | 0.47 GB | Q2/P2; QC/alignment only; publication linkage remains ambiguous |
| `PRJNA816041` | Eight paired-end ChIPmentation runs | 7.24 GB | Technical benchmark, outside the biological atlas |

For `PRJNA1173904`, the primary account identified during preparation was a
2026 preprint; a peer-reviewed gonadal-lncRNA study also reused the data. This
distinction is recorded in `publication_registry.tsv` and must remain explicit
in interpretation without changing the dataset's technical availability.

The so-called 108-library compendium combines `PRJNA236156` (60),
`PRJNA312093` (24), and `PRJNA602708` (24). It is not a ninth study: counting
it separately would duplicate FASTQs and evidence.

## Deposits not selected for execution

- `PRJEB23991`: eight paired-end libraries redundant with the
  miracidium/sporocyst subset of `PRJNA219440`. Biological sizes and deposited
  bases match, although ENA counts reads differently between deposits.
- `SRP007636`: four H3K9ac runs redeposited under `PRJNA384372` and
  `PRJNA219440`; the pairs are recorded in `duplicate_map.tsv`.
- `PRJEB3199`: two public runs without resolved mark, stage, control,
  replication, or publication. They do not meet the inclusion rule.

## Metadata reconstruction

The reconstruction does not invent groups. It maps public titles and aliases
to fields explicitly required by HelixForge: condition, biological replicate,
target, layout, control, and reference. Each library's original description
remains in `design.tsv` for audit.

The highest-impact decisions were:

- `PRJNA236156`: `Fa_1`, `Fb_2/Fb_4`, and `Fb_5/Fb_8` remain biological units;
  each IP points to its corresponding UB0.
- `PRJNA602708`: the generic input is shared by replicate number. “Beads
  without antibody” libraries are QC negative controls, not MACS3 controls.
- `PRJNA1173904`: A/A2 and B/B2 are normalized as replicates 1 and 2,
  respectively; each gonad/sexual-experiment group uses its public input.
- `PRJNA1090249`: only the eight *S. mansoni* runs are retained.
- `PRJNA219440`: two schistosomulum runs have unresolved targets in the
  deposit and retain that label. No run from this study advances to peak
  calling without a public control.

## HelixForge compatibility and frozen execution decisions

No structural pipeline change was required. Version 1.1.0 supports single-
and paired-end reads, Bowtie2, BAM processing, MACS3, FRiP, consensus/IDR,
annotation, and tracks. Because the API accepts one `peak_type` per invocation,
mixed-mark studies were split into `narrow` and `broad` run sets. This is an
execution configuration, not a pipeline fork.

Every IP is linked to its input/unbound control explicitly through
`control_id`; neither file order nor filenames are used to infer a control.
H3K4me3, H3K27ac, H3K9ac, H3K23ac, and H4K12ac follow the `narrow` path.
H3K27me3, H3K9me3, H3K36me3, H4K20me1, and HP1 follow the `broad` path.

In `v1.1.0`, independent invocations rebuild the Bowtie2 index and repeat
alignment; the ChIP-seq API cannot yet accept a prebuilt index or external BAMs
for the peak-calling path. Each run set therefore uses its own result directory
and is executed in series to prevent races and overwrites. This computational
duplication was recorded as an operational limitation, not changed here.

Frozen parameters and policies:

- Bowtie2 `--very-sensitive`;
- minimum MAPQ 30 and SAM exclusion flag 2308;
- duplicates retained in the analytical BAM and `MACS3 --keep-dup all`, per
  the validated HelixForge policy; duplication metrics remain mandatory;
- MACS3 3.0.4, `q=0.01`;
- initial effective size `391422299`, explicitly a proxy for total assembly
  length rather than an estimate of mappability;
- no blacklist: no validated, frozen *S. mansoni* blacklist was available.

## Limits on automatic inference

- Version 1.1.0 performs raw FastQC but has no native ChIP-seq trimming.
- IDR applies only to exactly two biological narrowPeak replicates.
- Broad marks and three-replicate sets use consensus supported by at least
  two replicates.
- `PRJNA219440` and `PRJNA384372` produce descriptive QC/BAM/tracks only;
  calling peaks without public inputs would be methodologically weak.
- Differential binding was not enabled for the first execution. Contrasts
  require prespecification after initial QC/peak review, without retrospective
  reclassification of P1.

## Preparation batch policy

The preparation batch limit was 150 GB of raw FASTQ. `PRJNA1173904` and
`PRJNA236156` therefore formed dedicated batches, while `PRJNA312093`,
`PRJNA602708`, `PRJNA1090249`, and `PRJNA219440` could be staged together
(104.06 GB). This did not authorize adding workdirs without checking space:
BAMs, indexes, and tracks must be counted separately, and each run set keeps
an independent directory. `PRJNA384372` remains an optional P2 batch.

## Project–atlas boundary

Results are reviewed study by study. Only peaks or consensus artifacts that
pass QC, have compatible reference identity, and carry complete provenance
may later enter the ChIP-seq atlas layer. Q2 BAMs, technical benchmarks, and
ambiguous deposits are not automatically included.
