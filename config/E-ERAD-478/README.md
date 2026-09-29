# E-ERAD-478

- Runs paired-end incluídos: **120**
- Amostras biológicas: **60**
- Amostras com consolidação de runs/lanes: **60**
- Runs ENA excluídos: **20**

A planilha corrigida S11/S12 do artigo é a autoridade para inclusão: 120 runs, 60 amostras e dois runs por amostra. O desenho ajusta `batch` e testa sexo, estado de pareamento e tempo por contrastes explícitos.

O `sample_id` repetido em `metadata.csv` é intencional: o HelixForge processa cada run e consolida os FASTQs no nível da amostra antes da quantificação. Os contrastes estão congelados em `de_spec.json`; não são inferidos dos nomes.
