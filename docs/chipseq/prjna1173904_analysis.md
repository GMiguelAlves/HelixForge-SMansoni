# PRJNA1173904: article-method processing of ovary H3K4me3

## Scope and source

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

## Processing policy

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

## Adaptations relative to the article

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

## Review and interpretation

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
