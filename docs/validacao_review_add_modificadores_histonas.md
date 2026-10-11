# Validação dos candidatos `REVIEW_ADD`

## Resultado

Os 11 candidatos originalmente classificados como `REVIEW_ADD` foram avaliados
contra o proteoma e a tabela de ortologia de *Schistosoma mansoni* SM_V10/WBPS19,
UniProtKB, perfis Pfam/HMMER e BLASTP contra proteínas humanas revisadas.

| Decisão | Genes | Consequência |
|---|---|---|
| Promover ao núcleo | `Smp_084140`, `Smp_127180` | passam para `CORE_ADD` |
| Manter exploratório | `Smp_185180`, `Smp_316920` | permanecem em `REVIEW_ADD` |
| Excluir do painel direto | `Smp_024380`, `Smp_059920`, `Smp_064250`, `Smp_128520`, `Smp_128530`, `Smp_204920`, `Smp_329500` | registrados em `EXCLUDED_REVIEW`; `Smp_064250` também permanece como contexto transcricional |

## Critério

Um candidato só foi promovido quando a evidência combinou domínio compatível,
arquitetura/ortologia funcional e uma função conhecida sobre histonas ou marcas de
histonas. A presença isolada de PHD ou Tudor não foi considerada suficiente, pois
esses domínios também ocorrem em proteínas que interagem com DNA, RNA ou outras
proteínas.

## Casos promovidos

### `Smp_084140` — eraser, sirtuína

- proteína de 321 aa;
- domínio SIR2/PF02146 completo e altamente significativo;
- anotação UniProt como sirtuína de classe I;
- melhor hit humano: SIRT2, com 46,1% de identidade, 89% de cobertura e
  `E = 5,46 × 10⁻97`;
- ortólogos de outros esquistossomos também são anotados como desacetilases
  dependentes de NAD.

Decisão: `CORE_ADD / ERASER / SIRTUIN_DEACETYLASE`.

### `Smp_127180` — reader, TAF3-PHD

- proteína de 1.543 aa com PHD/PF00628;
- ortólogo de *S. haematobium* anotado como subunidade 3 de TFIID;
- melhor hit humano: TAF3;
- o PHD de TAF3 reconhece diretamente H3K4me3, demonstrado estruturalmente e por
  ensaios de afinidade.

Decisão: `CORE_ADD / READER / TAF3_PHD_H3K4ME3_READER`.

## Casos mantidos para revisão

### `Smp_185180`

O domínio PHD está confirmado, mas não há ortologia funcional suficiente para
atribuir leitura de H3K4me3, H3K4me0 ou outra marca. Permanece exploratório.

### `Smp_316920`

O modelo de 161 aa contém dois PHDs e alinha com módulos PHD de KAT6A/KAT6B e
DPF. Entretanto, o gene codifica praticamente apenas o módulo tandem e não possui
ortologia de proteína inteira confiável. Deve permanecer exploratório até revisão do
modelo gênico ou validação experimental da especificidade de ligação.

## Casos excluídos do painel direto

- `Smp_024380`: o domínio RMT2 também ocorre em GAMT/metiltransferases de
  substratos não histônicos; o melhor hit humano foi GAMT, não uma PRMT de histona.
- `Smp_059920`: tem apenas 210 aa e não apresenta domínio SET/PF00856; não há
  suporte para classificá-lo como EZH/H3K27 methyltransferase.
- `Smp_064250`: contém PHD, mas a ortologia dominante é TCF20. Pode ser mantido
  como contexto de regulação transcricional, não como reader direto.
- `Smp_128520` e `Smp_128530`: não apresentaram PF02373/JmjC, motivo catalítico
  HxD ou ortólogo curado compatível com demetilase de histona.
- `Smp_204920` e `Smp_329500`: não apresentaram o perfil Tudor/PF00567 nem
  ortologia funcional que sustente leitura de histona; `Smp_329500` é um modelo de
  apenas 70 aa.

## Métodos e parâmetros

- referência: SM_V10, WormBase ParaSite WBPS19;
- proteoma: `schistosoma_mansoni.PRJEA36577.WBPS19.protein.fa.gz`;
- ortologia: `schistosoma_mansoni.PRJEA36577.WBPS19.orthologs.tsv.gz`;
- HMMER 3.4, perfis Pfam PF00567, PF00628, PF00856, PF02146, PF02373,
  PF12796, PF13771 e PF27704;
- BLASTP 2.12.0+ contra UniProtKB/Swiss-Prot humano, `E ≤ 1e-5`, até dez hits;
- os 11 loci possuem um único transcrito proteico no FASTA WBPS19 usado.

A evidência detalhada está em
`provenance/histone_modifier_review_validation/review_add_validation.tsv`.

## Fontes externas

- WormBase ParaSite WBPS19: https://parasite.wormbase.org/ftp.html
- UniProtKB REST: https://rest.uniprot.org/
- InterPro/Pfam: https://www.ebi.ac.uk/interpro/
- reconhecimento de H3K4me3 pelo PHD de TAF3:
  https://pubmed.ncbi.nlm.nih.gov/18682226/
