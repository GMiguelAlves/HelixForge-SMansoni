# PRJNA597909

- Amostras biológicas: 20
- Runs técnicos: 20
- Modo recomendado: `full` (QC → Salmon → tximport → DESeq2)
- Consolidação: o HelixForge recebe uma linha por run e usa o mesmo `sample_id` para aparar cada run e consolidá-los antes de Salmon.

## Grupos

- `adult_control`: 5 amostra(s)
- `adult_pzq`: 5 amostra(s)
- `juvenile_control`: 5 amostra(s)
- `juvenile_pzq`: 5 amostra(s)

## Arquivos

- `runs.tsv`: URLs, MD5, bytes e acesso ENA por run.
- `samples.tsv`: uma linha por amostra biológica.
- `metadata.csv`: contrato run-level consumido pelo HelixForge.
- `pipeline_config.sh` e `user_settings.sh`: adaptadores portáteis para o contrato atual.
- `de_spec.json`: desenho/contrastes versionados, quando estimáveis.
