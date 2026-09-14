# PRJNA602528

- Amostras biológicas: 10
- Runs técnicos: 10
- Modo recomendado: `import` (QC → Salmon → tximport; sem DE)
- Consolidação: o HelixForge recebe uma linha por run e usa o mesmo `sample_id` para aparar cada run e consolidá-los antes de Salmon.

## Grupos

- `time_0h`: 1 amostra(s)
- `time_0p25h`: 1 amostra(s)
- `time_12h`: 1 amostra(s)
- `time_1h`: 1 amostra(s)
- `time_24h`: 1 amostra(s)
- `time_3h`: 1 amostra(s)
- `time_48h`: 1 amostra(s)
- `time_6h`: 1 amostra(s)
- `time_96h`: 1 amostra(s)
- `time_9h`: 1 amostra(s)

O arquivo público tem uma biblioteca por tempo. Um modelo com dez níveis seria saturado e não teria replicação residual; por isso `de_spec.blocked.json` é deliberadamente não executável.

## Arquivos

- `metadata/PRJNA602528/runs.tsv`: URLs, MD5, bytes e acesso ENA por run.
- `metadata/PRJNA602528/samples.tsv`: uma linha por amostra biológica.
- `metadata/PRJNA602528/metadata.csv`: contrato run-level consumido pelo HelixForge.
- `pipeline_config.sh` e `user_settings.sh`: adaptadores portáteis para o contrato atual.
- `de_spec.blocked.json`: bloqueio explícito da inferência não estimável.
