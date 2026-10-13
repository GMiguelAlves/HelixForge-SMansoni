# Preparação do núcleo ChIP-seq de *Schistosoma mansoni*

## Escopo desta entrega

Esta preparação organiza os dados públicos antes da execução. Ela não altera o
HelixForge, não executa análises no servidor e não incorpora resultados no
atlas. O núcleo usa a referência já congelada `SM_V10/WBPS19` e o contrato
nativo ChIP-seq do HelixForge `v1.1.0` (`062001ceb46d8a25c2415823b97d0e310658f8c0`).

O inventário foi consultado novamente no ENA. Em 12 de outubro de 2026, WBPS19
continuava sendo a versão mais recente exposta no diretório oficial do WormBase
ParaSite. A referência local existente, seus checksums e sua validação são
reutilizados; não há segunda cópia do genoma no pacote ChIP-seq.

## Estudos preparados

| Estudo | Escopo | FASTQ bruto | Decisão |
|---|---:|---:|---|
| `PRJNA1173904` | 68 paired-end; oito marcas em gônadas | 130,77 GB | Q1/P1, peak calling com inputs por grupo |
| `PRJNA236156` | 60 single-end; quatro marcas, cercárias/adultos | 119,01 GB | Q1/P1, peak calling com UB0 por unidade biológica |
| `PRJNA312093` | 24 single-end; H3K4me3/H3K27me3 | 34,70 GB | Q1/P1, sexo e ciclo de vida, dois replicados |
| `PRJNA602708` | 24 paired-end; H3K4me3/H3K27me3 | 13,73 GB | Q1/P1, plasma/tratamento, dois replicados |
| `PRJNA1090249` | 8 single-end de *S. mansoni*; HP1 | 6,19 GB | Q1/P1 especializado; quatro runs de *Drosophila* excluídos |
| `PRJNA219440` | 29 runs mistos; conjunto histórico | 49,43 GB | Q2/P1, somente QC/alinhamento |
| `PRJNA384372` | 3 single-end H3K9ac | 0,47 GB | Q2/P2, somente QC/alinhamento e vínculo bibliográfico ainda ambíguo |
| `PRJNA816041` | 8 paired-end ChIPmentation | 7,24 GB | benchmark técnico fora do atlas biológico |

Para `PRJNA1173904`, o relato primário localizado é um preprint de 2026; há
também reutilização revisada por pares no trabalho de lncRNAs de gônadas. Essa
distinção permanece explícita no `publication_registry.tsv` e deverá constar da
interpretação final, sem alterar a disponibilidade técnica do conjunto.

O chamado compêndio de 108 bibliotecas é a agregação de `PRJNA236156` (60),
`PRJNA312093` (24) e `PRJNA602708` (24). Ele não é tratado como um nono estudo,
evitando duplicação de FASTQs e de evidência.

## Depósitos não executados

- `PRJEB23991`: oito bibliotecas paired-end redundantes com o bloco
  miracídio/esporocisto de `PRJNA219440`; os tamanhos biológicos e bases
  depositadas coincidem, embora o ENA conte reads de forma diferente entre os
  depósitos.
- `SRP007636`: quatro runs H3K9ac redepositados em `PRJNA384372` e
  `PRJNA219440`; os pares são registrados em `duplicate_map.tsv`.
- `PRJEB3199`: dois runs públicos, porém sem marca, estágio, controle,
  replicação ou publicação resolvidos. Não satisfaz a regra de inclusão.

## Reconstrução dos metadados

A reconstrução não inventa grupos. Ela transforma títulos e aliases públicos
em campos que o HelixForge exige explicitamente: condição, réplica biológica,
alvo, layout, controle e referência. O texto original de cada biblioteca é
mantido em `design.tsv` para auditoria.

As decisões de maior impacto são:

- `PRJNA236156`: `Fa_1`, `Fb_2/Fb_4` e `Fb_5/Fb_8` permanecem como unidades
  biológicas, e cada IP aponta para o UB0 correspondente.
- `PRJNA602708`: o input genérico é compartilhado por número de réplica; as
  bibliotecas “beads without antibody” ficam como controles negativos de QC,
  não como controle do MACS3.
- `PRJNA1173904`: A/A2 e B/B2 são normalizados, respectivamente, como réplicas
  1 e 2; cada grupo de gônada/experiência sexual usa seu input público.
- `PRJNA1090249`: somente os oito runs de *S. mansoni* são mantidos.
- `PRJNA219440`: dois runs de esquistossômulos têm alvo não resolvido no
  depósito e permanecem explicitamente rotulados como tal; nenhum run desse
  projeto avança para picos sem controle público.

## Compatibilidade com o HelixForge

Não é necessária mudança estrutural. A versão 1.1.0 aceita single- e paired-end,
Bowtie2, processamento de BAM, MACS3, FRiP, consenso/IDR, anotação e tracks.
Como a API atual recebe um único `peak_type` por invocação, estudos com marcas
misturadas foram divididos em entradas `narrow` e `broad`. Essa separação é uma
configuração de execução, não um fork do pipeline.

Em `v1.1.0`, cada invocação independente também reconstrói seu índice Bowtie2 e
refaz o alinhamento; a API ChIP-seq ainda não aceita um índice pré-construído ou
BAMs externos como entrada do caminho de picos. Para evitar corrida e
sobrescrita, cada run set usa diretório próprio e é executado em série. Essa
duplicação computacional é registrada como limitação operacional, sem mudar o
pipeline nesta entrega.

Parâmetros congelados:

- Bowtie2 `--very-sensitive`;
- MAPQ mínimo 30 e exclusão SAM 2308;
- duplicatas mantidas no BAM analítico e `MACS3 --keep-dup all`, conforme a
  política validada do HelixForge; métricas de duplicação continuam obrigatórias;
- MACS3 3.0.4, `q=0,01`;
- tamanho efetivo inicial `391422299`, documentado como proxy do comprimento
  total da montagem, não como estimativa de mapeabilidade;
- sem blacklist, porque ainda não existe blacklist de *S. mansoni* validada e
  congelada no projeto.

## Limitações que impedem inferência automática

- A versão 1.1.0 faz FastQC bruto, mas não implementa trimming ChIP-seq nativo.
- IDR só é usado para exatamente dois replicados biológicos narrowPeak.
- Marcas broad e conjuntos com três replicados usam consenso por suporte de
  pelo menos dois replicados.
- `PRJNA219440` e `PRJNA384372` produzem somente QC/BAM/tracks descritivos; sem
  inputs públicos, picos seriam metodologicamente frágeis.
- Differential binding não será acionado nesta primeira execução. Contrastes
  devem ser pré-especificados depois da revisão dos primeiros QC/peaks, sem
  reclassificação retrospectiva de P1.

## Política de lotes

O limite operacional desta entrega é 150 GB de FASTQ bruto por lote. Assim,
`PRJNA1173904` e `PRJNA236156` formam lotes dedicados, enquanto
`PRJNA312093`, `PRJNA602708`, `PRJNA1090249` e `PRJNA219440` podem ser
coestagiados (104,06 GB). Isso não autoriza somar workdirs sem verificar o
espaço: BAMs, índices e tracks são contabilizados separadamente e cada run set
mantém diretório independente. `PRJNA384372` permanece em um lote P2 opcional.

## Separação entre projeto e atlas

Os resultados serão revisados estudo a estudo. Somente picos/consensos que
passarem QC, tiverem identidade de referência compatível e provenance completo
serão promovidos posteriormente para a camada ChIP-seq do atlas. BAMs de Q2,
benchmarks técnicos e depósitos ambíguos não entram automaticamente.
