# PRJNA1173904: eight-mark gonad processing and annotation

## Full panel: accepted with warnings

The [full result panel](../../results/chipseq/PRJNA1173904/full_panel/README.md)
includes all **64 IP libraries** from
[Augusto, Quack, Grevelding and Grunau (2025)](https://doi.org/10.3390/ncrna11020025):
H3K4me3, H3K27me3, H3K9me3, H3K9ac, H3K23ac, H4K20me1, H4K12ac and H3K36me3.
Each mark has two labelled replicates in mating-experienced/inexperienced
ovaries and testes. Four Inputs were assessed separately. The reference is
SM_V10/WBPS19; this is an adapted processing reanalysis, not an exact historical
reproduction or the article's downstream combinatorial-state analysis.

### Processing

FastQC 0.12.1 precedes Bowtie2 2.5.5 end-to-end sensitive alignment with `gbar=4`.
The complete paired-end FASTQs are aligned without additional trimming or
downsampling in this full panel. Selection retains both primary, proper-pair,
QC-pass mates with MAPQ >=30, MAPQ !=255 and no `XS` tag. If either mate fails,
the pair is not retained. SAMtools 1.23.1 `fixmate -m` and `markdup -r` remove
duplicates per library.

ChromstaR 1.28.0 / R 4.3.3 fits an independent univariate HMM to each IP,
using 500 bp bins, 250 bp steps, paired-end counting and its native posterior
cutoff 0.5. Other settings include `short=TRUE`, standard initialization,
eps 0.1, at most 5,000 iterations, no chromosome prefit, count cap 500 and
seed 20261107. Upstream deduplication is not repeated inside ChromstaR.
The default maximum fragment width is 1,000 bp. Exact settings are in
[parameters.json](../../results/chipseq/PRJNA1173904/full_panel/manifests/parameters.json).
The posterior cutoff and BED score are not MACS p/q values.

Each native interval is retained separately. No shared-core intersection,
regional union, arbitrary gap filling or automatic replicate pooling changes
the regions. ChromstaR's exported `chr` aliases are explicitly mapped back to
the reference names without changing zero-based half-open BED coordinates.
Native exports, final per-peak posterior values and fit figures are preserved;
heavy fitted objects and bin-level arrays are not part of the retained full panel.

### Input policy

The 32 ovary IPs use their corresponding condition Input:
`SRR31021634` for experienced ovaries and `SRR31021585` for inexperienced
ovaries. Their final retained depths are 738,660 and 3,210,221 pairs, respectively;
the former remains a depth/quality limitation rather than an equivalent control.

The testis Inputs `SRR31021603` and `SRR31021623` have alignment rates of 7.25%
and 1.48% and retain only 206,169 and 55,731 pairs after filtering/deduplication.
A random screen of their unmapped reads found exogenous-sequence evidence;
that conditional small sample does not estimate whole-library contamination.
The current study deposit contains no additional testis Input.

These two Inputs are preserved but restricted from the primary background.
The **32 testis IPs are analyzed without Input**, explicitly as an alternative
branch. No ovary or unrelated-study control is substituted. Differences between
ovaries and testes must not be interpreted as controlled differential enrichment
under an identical background treatment.

### Annotation and interpretation

HOMER 4.11 annotates each library's intervals using the matching SM_V10 FASTA
and WBPS19 GTF with `-gid -gtf`. This is HOMER's native nearest-promoter
assignment, not proof of a regulated gene. Native classes/distances and the
source library are retained. The existing HF adapter transports the results
and checks nearest-TSS distance ties against the same HOMER TSS snapshot.
Its ambiguity fields do not assess every possible regulatory ambiguity.

The tables retain original/effective Input, mark, condition, replicate, caller,
native peak ID/score, library depth, raw QC flags, model genome fraction,
review warnings and `occupancy_status=not_assessed`. No per-library call or
gene assignment is silently upgraded to reproducible occupancy.

All 64 models reached their native convergence criterion and all 64 annotation
sets passed identity, coordinate, count and integrity checks. The 1,923,296
regions and nearest-promoter associations are sums across libraries, not unique
loci/genes. Numerical convergence and successful annotation do not establish
scientific validity of every interval.

### Included discordant libraries

| Group | Main warning |
|---|---|
| H3K4me3, inexperienced testis | SRR31021631 has 118,281 retained pairs and 48.50% genomic coverage versus 8,824,818 pairs and 6.36% in SRR31021630; basewise Jaccard is 0.103 |
| H3K4me3, experienced testis | Coverage differs between replicates (8.40% versus 3.39%); Jaccard is 0.350 |
| H3K9me3, experienced testis | Depth differs (1,595,924 versus 22,016,441 pairs), coverage is 14.31% versus 70.49%, and Jaccard is 0.197 |

The sparse SRR31021631 H3K4me3 fit has component means approximately 0.07 and
0.25 counts/bin; its broad inferred state is not independently established
enrichment. Its intervals remain available with a specific warning, not as
validation of its replicate. Other library QC warnings also remain in the tables.
No threshold was changed to increase concordance or match published counts.
Jaccard is descriptive overlap, not IDR, significance or a biological acceptance
test; broad genomic coverage can itself increase overlap.

**All libraries are accepted for inclusion with warnings, not excluded.**
Annotation acceptance does not authorize unqualified inferential comparisons.
Differential Binding remains **HOLD** pending appropriate control, replicate,
biological-unit and model review. No paired design is inferred from replicate
labels. The retained full panel is a results-only product: intervals, annotations,
figures, metrics, final peak posteriors and compact provenance. It does not retain
BAM/BAI, FASTQs, workdirs or fitted RDS objects. A future DB analysis would require
new alignments; the archive does not promise BAM-based re-entry.
No new DB, IDR, combinatorial-state or atlas result is claimed.

## Earlier H3K4me3 ovary calibration panel

The sections below describe only the earlier four-IP calibration, which remains
separate from the full panel and must not supply its control or depth policy.

### Scope and source

The accepted [public result panel](../../results/chipseq/PRJNA1173904/README.md)
follows the H3K4me3 ovary processing described by
[Augusto, Quack, Grevelding and Grunau (2025)](https://doi.org/10.3390/ncrna11020025),
*Chromatin Structure Around Long Non-Coding RNA (lncRNA) Genes in Schistosoma
mansoni Gonads*. It includes four IP libraries, two per mating experience,
and one common Input. The [sample table](../../results/chipseq/PRJNA1173904/metadata/sample_table.tsv)
identifies every run, control, sampling seed and selected-fragment count.

This result adopts the paper's processing approach, **not its downstream
eight-mark chromatin-state analysis**. It does not claim that every library
in the BioProject was reprocessed with ChromstaR or that the historical peak
catalogue was recovered. Earlier HF processing attempts are not part of the
accepted or retained result set and must not be relabelled as this
method-following panel.

### Processing policy

Existing QC and Bowtie2 2.5.5 alignments were reused against **SM_V10/WBPS19**.
Alignment selection retains MAPQ >=30. SAMtools 1.23.1 performs
template-coordinate duplicate removal before selection of complete fragments.
Each IP contributes exactly 2,500,000 randomly selected complete fragments,
keeping both mates; seeds are fixed. The common control `SRR31021585`
contributes all 4,509,411 eligible fragments.

| Tool | Fixed processing settings |
|---|---|
| MACS2 2.2.7.1 | BAM, one original R1 alignment per selected fragment; genome size 350,000,000; q=0.05; mfold 5–50; bandwidth 300; fragment model enabled |
| Peakranger 1.18 | CCAT low-resolution calling; extension 300 bp; FDR 0.05 |
| ChromstaR 1.28.0 / R 4.3.3 | Paired-end BAMs with Input; bins 500 bp; step 250 bp; mode `differential`; seed 2026110201; up to 5,000 iterations; univariate eps 0.1; multivariate eps 0.01 |

No second duplicate-removal pass was applied inside the callers after the
explicit upstream deduplication. MACS2 models succeeded without a `nomodel`
fallback; their estimated fragment lengths were 293, 307, 307 and 309 bp in
the sample-table IP order. Exact parameters are distributed in
[`parameters.json`](../../results/chipseq/PRJNA1173904/manifests/parameters.json).

The native ChromstaR posterior cutoff of 0.5 is not a MACS2 q-value or a
shared FDR threshold. Its configured mode does not make the interval-overlap
comparisons a Differential Binding analysis.

### Adaptations relative to the article

- SM_V10/WBPS19 replaces the article's v7 reference. Existing trimming,
  Bowtie2 version and MAPQ selection are retained; their equivalence to the
  historical processing is not established.
- The article's “2.5 million reads” is operationalized as **complete
  fragments**, not assumed to be an identical historical counting unit.
- The common single-sex ovary Input is consistent with the published control
  description; exact historical file identity is not independently proven.
  It is a control, not an additional biological replicate.
- Peakranger 1.18 replaces the cited 1.17. CCAT and its explicit FDR setting
  operationalize details not fully specified in the article.
- The MACS2 and ChromstaR versions are fixed here; the article does not fully
  specify their versions. The result is an adapted reanalysis, not an exact
  executable reproduction.

### Review and interpretation

All four univariate fits and the joint fit met native numerical convergence
criteria. The four native distribution plots are included with the results.
Replicate comparisons use only separately fitted univariate peaksets: joint
inference can make per-replicate exports identical and is not independent
evidence of replicate agreement.

Jaccard between experienced-ovary replicates is 0.712139 for ChromstaR,
0.685864 for MACS2 and 0.579643 for Peakranger. For inexperienced ovaries,
the values are 0.621837, 0.549628 and 0.555583, respectively. ChromstaR's
greater coverage contains 99.83–99.99% of MACS2-covered bases. Neither higher
overlap nor a larger catalogue demonstrates greater biological accuracy.

The published BED3 views preserve native interval boundaries. ChromstaR
sequence aliases are normalized through the explicit reference-name mapping;
coordinates remain zero-based and half-open. Union coverage is used only for
descriptive metrics, not to replace the per-library intervals with a consensus.

The panel does not generate gene associations or certify reproducible
occupancy. No paired statistical contrast is inferred from matching replicate
numbers, and no inferential DB or IDR result is published. The accepted scope
is the reviewed processing panel; incorporation into the atlas is separate.
