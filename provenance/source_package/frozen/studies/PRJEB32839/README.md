# PRJEB32839

- Amostras biológicas: 75
- Runs técnicos: 150
- Modo recomendado: `full` (QC → Salmon → tximport → DESeq2)
- Consolidação: o HelixForge recebe uma linha por run e usa o mesmo `sample_id` para aparar cada run e consolidá-los antes de Salmon.

## Grupos

- `cercaria_female`: 5 amostra(s)
- `cercaria_male`: 5 amostra(s)
- `cercaria_mixed`: 5 amostra(s)
- `egg_mixed`: 5 amostra(s)
- `juvenile_26d_female`: 5 amostra(s)
- `juvenile_26d_male`: 5 amostra(s)
- `miracidium_mixed`: 5 amostra(s)
- `schistosomulum_2d_female`: 5 amostra(s)
- `schistosomulum_2d_male`: 5 amostra(s)
- `schistosomulum_2d_mixed`: 5 amostra(s)
- `sporocyst_1d_mixed`: 5 amostra(s)
- `sporocyst_32d_female`: 5 amostra(s)
- `sporocyst_32d_male`: 5 amostra(s)
- `sporocyst_32d_mixed`: 5 amostra(s)
- `sporocyst_5d_mixed`: 5 amostra(s)

O ENA lista hoje 180 runs. Este pacote congela os 150 runs presentes no metadata oficial dos autores; os 30 restantes estão registrados em `excluded_ena_runs.tsv`.

## Arquivos

- `runs.tsv`: URLs, MD5, bytes e acesso ENA por run.
- `samples.tsv`: uma linha por amostra biológica.
- `metadata.csv`: contrato run-level consumido pelo HelixForge.
- `pipeline_config.sh` e `user_settings.sh`: adaptadores portáteis para o contrato atual.
- `de_spec.json`: desenho/contrastes versionados, quando estimáveis.
