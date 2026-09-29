# PRJEB3190

- Runs paired-end incluídos: **20**
- Amostras biológicas: **20**
- Amostras com consolidação de runs/lanes: **0**
- Runs ENA excluídos: **36**

Inclui somente a série poly(A)/mRNA canônica: 0 h tem 2 réplicas e os outros seis tempos têm 3 (20 amostras). Bibliotecas miRNA/small RNA e três rótulos não canônicos foram excluídos do caminho `full_length` e permanecem no inventário de exclusões.

Os valores `sequencing_batch_10363` e `sequencing_batch_10376` foram derivados dos prefixos de `experiment_alias` da ENA. A ENA define alias como um nome fornecido pelo submitter, e o artigo não documenta esses prefixos como lotes técnicos. Por isso eles são preservados para PCA, diagnóstico e eventual investigação, mas não são usados como covariável inferencial. O desenho congelado permanece `~ condition`, acompanhado desse aviso.

O `sample_id` repetido em `metadata.csv` é intencional: o HelixForge processa cada run e consolida os FASTQs no nível da amostra antes da quantificação. Os contrastes estão congelados em `de_spec.json`; não são inferidos dos nomes.
