# Download e referência

## Referência

Submeta `scripts/download_reference.sh` em um nó de compute com rede, ou execute-o dentro de uma alocação interativa curta. Ele baixa genome, mRNA transcripts, GFF3 e GTF canônico WBPS19, verifica gzip, descompacta, calcula SHA-256 e grava `REFERENCE_READY.json`.

## FASTQs

Para cada projeto:

```bash
bash scripts/submit_downloads.sh PRJEB32839
bash scripts/submit_downloads.sh PRJEB14695
bash scripts/submit_downloads.sh PRJNA597909
bash scripts/submit_downloads.sh PRJNA602528
```

Cada tarefa baixa as duas mates, valida bytes, MD5 e gzip e grava um marcador por run. Valide novamente após os arrays:

```bash
python3 scripts/validate_downloads.py PRJEB32839 --package-root "$HF_PACKAGE_ROOT" --data-root "$HF_DATA_ROOT"
```

Repita para os quatro projetos. Não concatene FASTQs manualmente: o contrato atual recebe cada run em uma linha e o módulo de QC consolida os runs que compartilham `sample_id` depois do trimming.
