# HelixForge: quatro estudos RNA-seq de Schistosoma mansoni

Pacote operacional para PRJEB32839, PRJEB14695, PRJNA597909 e PRJNA602528, alinhado ao HelixForge `8a1eac201d34861a9f8d3cd595b472a5ca76cc97` e à referência `Schistosoma_mansoni_SM_V10_WBPS19`.

## Resultado da revisão

Não foi necessária alteração estrutural no HelixForge. Os inputs seguem o contrato run-level existente; runs técnicos compartilham `sample_id` e são consolidados pelo QC antes de Salmon. Salmon, tximport (`lengthScaledTPM`) e DESeq2 permanecem no caminho certificado.

WBPS19 (release 19, março de 2024) é a versão mais recente exibida pelo WormBase ParaSite para SM_V10 na consulta de 11/09/2026. O pacote fixa genome, mRNA transcripts, GFF3 e GTF canônico dessa mesma release.

## Escopo congelado

| Projeto | Runs | Amostras | Execução |
|---|---:|---:|---|
| PRJEB32839 | 150 | 75 | full |
| PRJEB14695 | 138 | 23 | full |
| PRJNA597909 | 20 | 20 | full |
| PRJNA602528 | 10 | 10 | import, sem DE |

Total: 318 runs, 128 amostras e cerca de 501,1 GB de FASTQ comprimido. Reserve inicialmente cerca de 6 TB de scratch, coerente com o teto operacional de 10,56× usado no catálogo.

## Início rápido

1. Leia `runbooks/01_PREPARACAO_SERVIDOR.md`.
2. Configure e carregue `env/server.env`.
3. Baixe/valide a referência e os FASTQs.
4. Rode `scripts/preflight_server.sh`.
5. Comece pelo piloto `bash scripts/run_study.sh PRJNA597909`.

## Limitações preservadas

- PRJNA602528 não tem replicação pública por tempo; DESeq2 é bloqueado.
- A API DE 1.0 não expressa interação formal maturidade × tratamento; os quatro grupos de PRJNA597909 recebem contrastes pareados simples.
- Apptainer/Conda em Slurm dependem da configuração do site e permanecem experimentais na documentação pública do HelixForge; o piloto é obrigatório.
- O checkout local do usuário não foi alterado. Este pacote foi auditado contra uma cópia limpa da `master` pública atual.
