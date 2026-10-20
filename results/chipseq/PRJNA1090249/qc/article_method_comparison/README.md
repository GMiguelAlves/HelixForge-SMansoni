# HP1 — comparação descritiva de métodos

Abra [o relatório com figuras](relatorio.html). As três figuras estão em
PNG e PDF no diretório `figures/`. As tabelas preservam os valores usados.

## Conclusão

O processamento e a avaliação terminaram sem erros. A adaptação do método do
artigo à SM_V10/WBPS19 permitiu ajustar modelos após agrupar bibliotecas, mas
não estabeleceu um conjunto robusto de domínios de HP1. Os resultados permanecem
exploratórios, com avisos; Differential Binding não foi aprovado.

O aumento de cobertura próximo ao final dos genes também aparece no Input.
Essa forma do perfil, sozinha, não demonstra enriquecimento específico de HP1.

## Diferenças importantes

- SM_V10/WBPS19, não a referência antiga do artigo.
- Três bibliotecas IP por estágio; o subconjunto histórico exato não foi confirmado.
- Profundidade realizada: 3.819.229 reads por biblioteca, não 4,7 milhões.
- MACS2 2.2.9.1 e Picard 3.4.0; versões originais não estabelecidas.
- Parâmetros geométricos de broad calling explicitados; não ajustados para imitar contagens publicadas.
- Pools servem à descrição do sinal, não são novas réplicas biológicas.

Nenhum alinhamento, BAM original ou resultado científico principal foi alterado.
Este pacote não contém caminhos privados de execução ou logs operacionais.

Método e resultados resumidos também estão na
[documentação do estudo](https://github.com/GMiguelAlves/HelixForge-SMansoni/blob/analysis/prjna1090249-chipseq/docs/chipseq/prjna1090249_peak_policy.md).
