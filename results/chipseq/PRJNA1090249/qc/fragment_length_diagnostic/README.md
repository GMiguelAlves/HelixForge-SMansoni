# Fragment-length diagnostic — PRJNA1090249

Supplementary QC figures for all eight libraries: six HP1 IPs and two Inputs.
Open **index.html** for the figure gallery. Each library includes PNG/PDF plots
and the underlying strand-correlation TSV. The three available native MACS3
models are included as separate PDFs.

## Scope and interpretation

These are exploratory diagnostic outputs, not a completed ChIP-seq terminal
report, a certification of biological replication, or an estimate approved for
`--extsize`. Fragment-length assessment remains **INCONCLUSIVE**.

The observed modal read length is 76 bp. Native MACS3 models selected 71–72 bp
in the three cercarial libraries. A read-length-adjacent component is also
present in Inputs and is sensitive to exact-start deduplication. This raises a
possible technical-artifact concern but does not prove its cause, establish a
physical fragment length, or imply that the libraries are unusable.

No new fragment length, duplicate policy or peak-calling parameter was adopted
to produce these figures. BAMs were read without modification.

The subsequent [duplicate-policy model sensitivity review](model_sensitivity_review.md)
did not recover an estimable fragment model by limiting tags to one per
position/strand. This does not establish that duplicate handling explains the
difference from a historical successful run.

## Method and figure key

Genome-wide strand Pearson correlation in 5-bp bins, over displacements from
−500 to +1500 bp. All reference contigs, including mitochondrial sequence,
were included without a mappability/repeat mask. Unmapped, secondary,
supplementary and QC-fail alignments and MAPQ<30 were excluded.

- Solid blue: all eligible alignments.
- Dashed orange: unique 5-prime coordinates per strand, collapsed before binning.
- Gray band: modal read length ±10 bp.
- Dotted line: selected native MACS3 d, when available.
- Upper panel: full displacement range; lower panel: 0–500-bp zoom.

Each library has its own Y-axis scale. The curves are descriptive, without
inferential confidence intervals. Exact-start collapse is a diagnostic
sensitivity analysis; coincident alignments are not all proven PCR duplicates.
No official SPP/phantompeak NSC or RSC quality score was calculated.

The figures should accompany the final study results even if original
BioAnalyzer profiles are recovered or a different modeling policy is approved.
Keep their original diagnostic scope and do not relabel them as peak-calling
results from a subsequent parameterization.

## Decisão posterior: extensão fixa

A [política de peak calling](peak_policy.md) registra a extensão fixa de 200 bp, sem modelo automático. É uma hipótese operacional, não uma estimativa experimental. O processamento principal e as 12 chamadas de sensibilidade a 150 e 250 bp concluíram sem erros. A comparação demonstrou sensibilidade relevante à extensão; consulte [a tabela completa](extension_sensitivity.tsv). O estudo foi encerrado e aceito para inclusão no atlas, mantendo essas ressalvas quantitativas. DB não será executado nesta análise. Os diagnósticos anteriores não certificam um tamanho físico de fragmento nem fronteiras robustas dos domínios.
