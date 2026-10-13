# Metadados RNA-seq

Cada diretório de estudo contém `runs.tsv` (aquisição e integridade dos FASTQs),
`samples.tsv` (unidades biológicas) e `metadata.csv` (entrada do HelixForge).
Exclusões explícitas do depósito público, quando necessárias, ficam em
`excluded_ena_runs.tsv`.

`master_sample_manifest.tsv` é o índice conservador dos estudos inicialmente
consolidados. `p1_pending_manifest.tsv` registra um lote histórico de preparação;
não define um projeto científico separado nem substitui os registros individuais
dos estudos. Os arquivos em `vocabularies/` preservam os valores curados sem
presumir equivalência biológica entre estudos.

Os contratos científicos e o estado de cada estudo estão descritos no
[registro de estudos](../../docs/studies_and_preparation.md).
