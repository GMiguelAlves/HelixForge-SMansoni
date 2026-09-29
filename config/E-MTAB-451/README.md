# E-MTAB-451

- Runs paired-end incluídos: **12**
- Amostras biológicas: **11**
- Amostras com consolidação de runs/lanes: **1**
- Runs ENA excluídos: **0**

DE restrita aos três grupos replicados. `adult_mixed` e `cercarial_tail` (n=1) permanecem para quantificação/exploração, mas não aparecem em contrastes inferenciais. `somule1` possui dois runs técnicos e será consolidada.

O `sample_id` repetido em `metadata.csv` é intencional: o HelixForge processa cada run e consolida os FASTQs no nível da amostra antes da quantificação. Os contrastes estão congelados em `de_spec.json`; não são inferidos dos nomes.
