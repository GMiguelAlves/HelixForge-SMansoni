# PRJEB14695

- Amostras biológicas: 23
- Runs técnicos: 138
- Modo recomendado: `full` (QC → Salmon → tximport → DESeq2)
- Consolidação: o HelixForge recebe uma linha por run e usa o mesmo `sample_id` para aparar cada run e consolidá-los antes de Salmon.
- Concorrência: no máximo 8 jobs simultâneos no Slurm, tanto na aquisição quanto no workflow.
- Unidade estatística: DESeq2 recebe 23 amostras biológicas; os 138 runs técnicos nunca são tratados como réplicas independentes.

## Grupos

- `ovary_mixed_sex`: 3 amostra(s)
- `ovary_single_sex`: 2 amostra(s)
- `testis_mixed_sex`: 3 amostra(s)
- `testis_single_sex`: 3 amostra(s)
- `whole_female_mixed_sex`: 3 amostra(s)
- `whole_female_single_sex`: 3 amostra(s)
- `whole_male_mixed_sex`: 3 amostra(s)
- `whole_male_single_sex`: 3 amostra(s)

## Arquivos

- `metadata/PRJEB14695/runs.tsv`: URLs, MD5, bytes e acesso ENA por run.
- `metadata/PRJEB14695/samples.tsv`: uma linha por amostra biológica.
- `metadata/PRJEB14695/metadata.csv`: contrato run-level consumido pelo HelixForge.
- `pipeline_config.sh` e `user_settings.sh`: adaptadores portáteis para o contrato atual.
- `de_spec.json`: desenho/contrastes versionados, quando estimáveis.
