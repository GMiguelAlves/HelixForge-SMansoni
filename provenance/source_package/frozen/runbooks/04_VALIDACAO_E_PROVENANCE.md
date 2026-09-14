# Validação e provenance

Critérios mínimos antes de interpretar biologia:

- 318 runs paired-end selecionados e 128 amostras biológicas no total.
- PRJEB32839 usa apenas as 75 amostras publicadas (150 runs). Os 30 runs adicionais atualmente associados ao projeto ENA ficam excluídos e documentados.
- PRJEB14695 mantém 138 runs em 23 amostras, seis runs técnicos por amostra.
- PRJNA597909 mantém 20 amostras independentes, cinco por grupo.
- PRJNA602528 tem uma biblioteca por tempo e não recebe inferência DESeq2.
- todos os downloads passam bytes + MD5 + gzip; referência passa SHA-256 e verificações estruturais.
- versões/commit, comandos, configs, trace e manifests terminais são preservados.

Ao final, execute `bash scripts/collect_provenance.sh`. Não mova `rnaseq_run_manifest.json` sem o diretório `integration_artifacts/` correspondente.
