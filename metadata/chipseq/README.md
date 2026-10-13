# Metadados ChIP-seq

Este diretório registra o núcleo público de ChIP-seq de *Schistosoma mansoni*
sem misturá-lo aos estudos RNA-seq já existentes. O registro não promove
automaticamente nenhum resultado ao atlas.

## Arquivos globais

- `study_registry.tsv`: decisão de inclusão, prioridade, escopo executável e
  volume bruto por depósito.
- `run_inventory.tsv`: inventário consolidado de FASTQs públicos.
- `publication_registry.tsv`: vínculo entre depósitos e publicações, incluindo
  os casos ambíguos.
- `duplicate_registry.tsv` e `duplicate_map.tsv`: depósitos que não devem ser
  baixados novamente.

## Arquivos por estudo

- `runs.tsv`: URLs HTTPS do ENA, MD5, tamanho e contagens para aquisição
  verificável.
- `design.tsv`: reconstrução curada de estágio, sexo, tratamento, unidade
  biológica, réplica, alvo e controle.
- `helixforge.narrow.tsv` e `helixforge.broad.tsv`: entradas nativas do
  HelixForge separadas por geometria de pico.
- `helixforge.alignment.tsv`: entrada restrita a QC/alinhamento quando não há
  controle público adequado.

Os FASTQs não são versionados. Seus nomes locais em `helixforge.*.tsv`
correspondem ao nome final de cada URL em `runs.tsv`.

## Regras científicas congeladas

1. Inputs/unbound são ligados explicitamente por `control_id`; controles não
   são inferidos por ordem ou nome de arquivo.
2. H3K4me3, H3K27ac, H3K9ac, H3K23ac e H4K12ac usam o caminho `narrow`.
3. H3K27me3, H3K9me3, H3K36me3, H4K20me1 e HP1 usam o caminho `broad`.
4. `PRJNA219440` e `PRJNA384372` não avançam para peak calling, pois os
   controles públicos necessários não estão disponíveis.
5. `PRJNA816041` é um benchmark técnico de ChIPmentation e não integra a
   evidência biológica do atlas.
6. `PRJEB23991` e `SRP007636` são redundantes com depósitos canônicos e não
   devem ser baixados.

As consultas ENA completas que sustentam esses arquivos estão congeladas em
`provenance/source_package/frozen/chipseq/<accession>/ena_read_run.tsv`.
