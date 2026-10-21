# HelixForge-SMansoni RNA-seq Atlas

Open `atlas.html` for RNA-seq or `chipseq.html` for ChIP-seq. Both explorers share gene collections, including the 118-gene machinery preset. `machinery.html` redirects to that collection. The atlas integrates 12 accepted studies, 268 biological samples, 535 technical runs, 9,914 shared genes, and 71 original per-study DESeq2 contrasts.

The reviewed HF ChIP-seq layer includes 2 study, 22 deposited IP libraries, 244788 regional candidates and 9523 genes with native associations. Gene associations and loci are available throughout the atlas gene universe. HOMER nearest-promoter assignments include source-classified intergenic regions; location class and nearest-gene assignment remain distinct.

The gene explorer includes an SM_V10 coordinate-aware locus view of every native gene–segment association. The display is descriptive and does not imply sample pairing or causality.

Separate literature layers contain 9,516 Anderson annotation/source claims and 154 Wang class-dependent sporocyst stem-cell transcripts. These are published evidence, not extra bulk samples, atlas TPM values or DESeq2 contrasts. All source records, including unmapped IDs, are available in the external evidence TSVs.

The literature annotation layer additionally contains 7,431 Vasconcelos lncRNA models (181 robust candidates; 916 with published stage-DE evidence) and 11,907 Almeida historical transcript claims. Vasconcelos coordinates remain labelled as v5.2 source coordinates, and Almeida contigs without a stable gene ID remain unmapped.

The gene explorer shows the current SM_V10/WBPS19 annotation for every atlas gene. Historical model claims cannot be confirmed from a stable-ID match alone.

The atlas is descriptive. It does not fit a joint differential-expression model, perform batch correction, or reinterpret the original contrast orientations. Machine-readable tables are under `data/`, exportable figures under `figures/`, and all generated artifacts are recorded in `manifest.json`.
