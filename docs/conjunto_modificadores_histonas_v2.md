# Conjunto de modificadores de histonas v2 (SM_V10/WBPS19)

## Decisão de curadoria

A lista histórica é uma curadoria biológica do projeto. Os membros de `Writers`,
`Readers` e `Erasers` são tratados como modificadores de histonas confirmados pela
curadoria da autora. O grupo HP1 é preservado, mas separado como leitor/organizador
de heterocromatina e não entra na contagem de modificadores enzimáticos.

O estado do identificador é uma questão diferente da classificação biológica. Os cinco
IDs históricos que antes estavam ambíguos ou ausentes foram resolvidos por comparação
de sequência e, nos loci fragmentados, pela preservação do segmento funcional esperado.

## Camadas da lista

| Camada | Critério | Genes únicos | Uso no atlas |
|---|---|---:|---|
| `CORE_CURATED` | lista histórica curada pelo projeto, excluindo HP1 da contagem enzimática | 60 | núcleo, após resolver IDs pendentes |
| `CORE_ADD` | anotação explícita, catálogo publicado ou validação de domínio+ortologia | 38 | expansão principal |
| `REVIEW_ADD` | domínio confirmado, mas função sobre histona ou modelo gênico ainda inconclusivo | 2 | painel exploratório, fora do núcleo |
| `CONTEXT` | HP1, cofatores, remodeladores e chaperonas | 24 | contexto mecanístico, fora da contagem de modificadores diretos |
| `EXCLUDED_REVIEW` | hipótese inicial não sustentada como modificador/reader direto | 7 | documentado, mas não analisado como candidato |

Há **126 genes de entrada únicos** em todas as camadas. A tabela é
normalizada por associação funcional: um gene pode aparecer em mais de uma linha,
como uma acetiltransferase que também contém bromodomínio.

## Classes importantes acrescentadas

- **Writers de acetilação:** HATs explícitas em SM_V10, ausentes da lista histórica.
- **Erasers de acetilação:** HDACs clássicas e sirtuínas; `Smp_084140` foi
  promovido após confirmação do domínio SIR2 e ortologia com SIRT2.
- **Readers de acetilação:** 22 proteínas do catálogo publicado de bromodomínios.
- **Reader de H3K4me3:** `Smp_127180` foi promovido como ortólogo de TAF3 com
  PHD; a função de leitura de H3K4me3 do PHD de TAF3 é experimentalmente estabelecida.
- **PHDs ainda exploratórios:** `Smp_185180` e `Smp_316920` têm domínio PHD,
  mas ainda não possuem suporte suficiente para especificidade de marca de histona.
- **Hipóteses rejeitadas para o painel direto:** `Smp_024380`, `Smp_059920`,
  `Smp_064250`, `Smp_128520`, `Smp_128530`, `Smp_204920` e `Smp_329500` foram
  mantidos como registros de auditoria em `EXCLUDED_REVIEW`.
- **Contexto cromatínico:** Polycomb/Trithorax, CoREST/SAP18/RBBP4, remodeladores
  ATP-dependentes e chaperonas de histonas são mantidos em grupos funcionais próprios.

## Contagens por camada e classe

| Camada | Classe | Associações gene-classe |
|---|---|---:|
| `CORE_CURATED` | `ERASER` | 13 |
| `CORE_CURATED` | `READER` | 24 |
| `CORE_CURATED` | `WRITER` | 26 |
| `CORE_ADD` | `ERASER` | 9 |
| `CORE_ADD` | `READER` | 23 |
| `CORE_ADD` | `WRITER` | 10 |
| `REVIEW_ADD` | `READER` | 2 |
| `CONTEXT` | `CHAPERONE` | 4 |
| `CONTEXT` | `COFACTOR` | 10 |
| `CONTEXT` | `HP1` | 2 |
| `CONTEXT` | `REMODELER` | 8 |
| `EXCLUDED_REVIEW` | `ERASER` | 2 |
| `EXCLUDED_REVIEW` | `READER` | 3 |
| `EXCLUDED_REVIEW` | `WRITER` | 2 |


## Estado dos IDs da curadoria histórica

- correspondência exata SM_V10: 48 associações;
- correspondência única por ID anterior: 12 associações;
- correspondência exata por sequência: 2 associações;
- correspondência de alta confiança por sequência: 1 associações;
- sucessor funcional de locus fragmentado: 2 associações;
- correspondência ambígua: 0 associações;
- ainda não encontrado em WBPS19: 0 associações.

Esses números são contagens de associações funcionais, não de evidência biológica.
A evidência de resolução dos cinco IDs está congelada em
`provenance/histone_modifier_id_resolution/historical_id_resolution.tsv`.

## O que não foi misturado ao núcleo

- **HP1:** reader/organizador de heterocromatina, mantido como contexto por decisão do projeto.
- **Histonas canônicas e variantes:** são substratos do sistema, não writers/erasers/readers.
- **Remodeladores e chaperonas:** alteram acesso/organização da cromatina, mas não
  necessariamente escrevem, apagam ou leem uma marca covalente.
- **Genes inferidos apenas por domínio:** somente os dois PHDs inconclusivos
  permanecem em `REVIEW_ADD`; hipóteses contraditas pelos perfis/ortologia foram excluídas.

## Arquivos e ativação

A tabela completa e congelada está em `config/candidates/histone_modifier_gene_set_v2.tsv`.
O arquivo histórico `config/PRJEB32839/candidate_genes.txt` não foi alterado, para
preservar a reprodutibilidade do atlas já construído. A v2 só deve substituir a lista
ativa somente nas análises que declararem explicitamente esta versão. Os cinco IDs
pendentes já foram resolvidos; os dois PHDs de `REVIEW_ADD` continuam exploratórios.

A auditoria dos antigos `REVIEW_ADD`, incluindo HMMER, BLASTP e ortologia, está em
`docs/validacao_review_add_modificadores_histonas.md` e
`provenance/histone_modifier_review_validation/review_add_validation.tsv`.

## Fontes

- Curadoria histórica: lista fornecida e confirmada pela autora do projeto.
- Referência publicada de contexto para a maquinaria de metilação: https://pmc.ncbi.nlm.nih.gov/articles/PMC6125318/
- Catálogo publicado de proteínas com bromodomínio: https://pmc.ncbi.nlm.nih.gov/articles/PMC10726355/
- Função leitora de H3K4me3 do PHD de TAF3: https://pubmed.ncbi.nlm.nih.gov/18682226/
- Anotação local: `config/reference/gene_annotations_wbps19.tsv` (SM_V10/WBPS19).
