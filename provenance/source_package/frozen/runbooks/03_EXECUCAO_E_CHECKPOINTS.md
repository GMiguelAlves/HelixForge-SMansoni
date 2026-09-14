# Execução e checkpoints

## Ordem recomendada

1. `bash scripts/preflight_server.sh`
2. PRJNA597909 como piloto pequeno.
3. PRJEB14695.
4. PRJEB32839.
5. PRJNA602528 em modo `import`.

Comando:

```bash
bash scripts/run_study.sh PRJNA597909
```

## Checkpoints por estudo

- Metadata: `pipeline_info/native_rnaseq/metadata/metadata_validation.json` deve indicar `status=valid`, o total de runs e amostras esperado.
- QC: FastQC/Trim Galore completos; o plano deve mostrar um `merged_sample_r1/r2` comum a todos os runs da amostra.
- Salmon: um `quant.sf` por amostra biológica e taxa de mapeamento revisada; investigar amostras muito abaixo da distribuição do estudo.
- Import: matrizes gene-level sem IDs ausentes; política `full_length + lengthScaledTPM`.
- DE (três estudos): preflight sem desenho rank-deficient e contrastes com pelo menos duas réplicas por nível.
- Terminal: para execuções `full`, localizar `rnaseq_run_manifest.json`, validar checksums e preservar o diretório irmão `integration_artifacts/`.
- Nextflow: revisar `trace.tsv`, `execution_report.html`, `timeline.html`, `.nextflow.log`, retries e jobs com exit code não zero.

Sempre retome com `-resume` e o mesmo `HF_WORK_ROOT`. O primeiro estudo constrói o índice Salmon; execuções seguintes podem reutilizar a tarefa pelo cache se referência, parâmetros e workdir forem idênticos.
