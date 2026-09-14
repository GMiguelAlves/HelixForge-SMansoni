# Auditoria do contrato HelixForge

Base revisada: branch pública `master`, commit `8a1eac201d34861a9f8d3cd595b472a5ca76cc97`.

- Entrada RNA-seq: metadata run-level com `dataset`, `sample_id` e `run_accession`; FASTQs paired são resolvidos em `<SCRATCH_ROOT>/<dataset>/fastq_ftp/<run>_{1,2}.fastq.gz`.
- Consolidação: runs que compartilham `sample_id` são aparados individualmente e reunidos nos artefatos `merged_sample_r1/r2` antes da quantificação.
- Referência Salmon: transcriptoma e GTF/GFF3 obrigatórios; genome também é fixado no pacote para identidade SM_V10 e uso futuro.
- Import: caminho production_v1 exige bibliotecas full-length com `countsFromAbundance=lengthScaledTPM`.
- DE: API 1.0 exige DESeq2/Wald, fórmula na ordem covariáveis + variável, desenho full-rank, contrastes explícitos e replicação mínima.
- Terminal: `rnaseq_run_manifest.json` é produzido no modo `full`; PRJNA602528 para em `import` porque DE não é estimável e, portanto, não produz o terminal completo.
- Slurm: os processos usam o executor Slurm; aquisição de dados permanece fora do DAG científico.

Conclusão: nenhuma mudança no core é necessária para os três estudos com DE estimável e para quantificação/import do quarto. As limitações são de desenho experimental e suporte estatístico, não de organização de arquivos.
