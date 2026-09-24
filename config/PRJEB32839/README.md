# PRJEB32839

- Amostras biológicas: 75
- Runs técnicos: 150
- Modo recomendado: `full` (QC → Salmon → tximport → DESeq2)
- Consolidação: o HelixForge recebe uma linha por run e usa o mesmo `sample_id` para aparar cada run e consolidá-los antes de Salmon.
- Runtime fixado: HelixForge `v1.0.2` (`5d4b3e696319db5cd7633472504964f1dc7c0434`), sem atualização automática para `master`.

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

- `metadata/PRJEB32839/runs.tsv`: URLs, MD5, bytes e acesso ENA por run.
- `metadata/PRJEB32839/samples.tsv`: uma linha por amostra biológica.
- `metadata/PRJEB32839/metadata.csv`: contrato run-level consumido pelo HelixForge.
- `pipeline_config.sh` e `user_settings.sh`: adaptadores portáteis para o contrato atual.
- `de_spec.json`: desenho/contrastes versionados, quando estimáveis.

## Gate de armazenamento

Nenhum FASTQ pode ser adquirido antes de `provenance/PRJEB32839/storage_preflight.json`
registrar `authorization: GO`. O cálculo usa o inventário congelado de 150 runs,
a soma exata dos 300 FASTQs e uma projeção conservadora de 6× com 20% de margem.
