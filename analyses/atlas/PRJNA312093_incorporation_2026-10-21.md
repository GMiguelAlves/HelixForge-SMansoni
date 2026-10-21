# Incorporação de ChIP-seq — 21/10/2026

Os catálogos publicados de PRJNA312093 foram incorporados ao atlas local, ao navegador de regiões, às páginas de genes e à integração de coleções com RNA-seq. Fonte congelada: HelixForge-SMansoni, commit `e583c84ef5cd37edd2f4df8450f17cc42c39c279`. A fonte HP1 anterior permanece preservada em `c9e77fee805d73e0527a02f83217c7a2763d270c`.

## Conteúdo incorporado

| Produto | Quantidade |
|---|---:|
| H3K27me3, catálogo regional broad | 188.268 regiões |
| H3K4me3, catálogo regional narrow | 53.837 regiões |
| Total novo | 242.105 regiões |
| Registros de suporte por biblioteca, incluindo suporte zero | 484.210 |
| Intervalos suplementares H3K4me3 com IDR ≤ 0,05 | 18.179 |
| Segmentos suplementares de núcleo compartilhado, novas marcas | 90.485 |
| Genes com atribuição na fonte PRJNA312093 | 9.504 |
| Atlas completo, incluindo HP1 | 244.788 regiões / 9.523 genes com atribuições |

Cada marca possui quatro contextos: cercárias e adultos, masculinos e femininos, com duas bibliotecas IP por contexto. São 16 IPs e oito Inputs únicos compartilhados entre as marcas. Incluindo HP1, há 22 IPs e dez Inputs únicos. O campo histórico `chipseq_controls=18` soma participações dos controles nos contextos, não bibliotecas únicas.

As coordenadas SM_V10/WBPS19, classificações, atribuições escolhidas pela fonte, IDs alternativos e estados de ambiguidade foram preservados. H3K4me3 contém 22 atribuições ambíguas, 52.644 sem empate de distância detectado e 1.171 não avaliadas. H3K27me3 mantém 188.268 registros não avaliados. Os intervalos IDR e núcleos compartilhados são camadas opcionais separadas do catálogo mínimo de uma biblioteca.

## Limites científicos preservados

- A admissão é descritiva, conforme a revisão da fonte. Differential Binding continua em HOLD e não foi executado.
- A comparação estrita com a ferramenta original permanece FAIL para uma atribuição em cada contexto masculino de H3K27me3, envolvendo genes com TSS compartilhado, Smp_095580 e Smp_341610. Não houve alteração para forçar concordância.
- Pools femininos selecionados permanecem sem resolução; sexo e instrumento estão completamente confundidos nas comparações de adultos. As demais ressalvas de QC estão na página do estudo.
- O pacote publicado não contém coordenadas de chamadas individuais de PRJNA312093 nem bigWigs. Essas trilhas não foram inventadas; a cobertura por biblioteca pode ser consultada nos detalhes das regiões.
- A integração usa estágio e sexo registrado quando aplicável; RNA misto não é equiparado a ChIP de sexo específico. Não há pareamento de amostras ou inferência causal.
- As atribuições HOMER ao promotor mais próximo, inclusive intergênicas, não estabelecem regulação ou ausência biológica de uma marca.
- Quatro genes associados não estão no universo RNA: Smp_334270, Smp_334590, Smp_334690 e Smp_335440, totalizando 23 regiões. Permanecem no ChIP e nos downloads; nenhuma expressão foi imputada.
- Permanecem exatamente os mesmos dados de RNA, metadados e proveniência das análises de RNA: 9.914 genes, 268 amostras e 71 contrastes. Os contrastes de adultos/cercárias de PRJNA312093 continuam excluídos.

## Entrega e tamanho

O carregamento inicial usa resumos compactos. Coordenadas dos novos catálogos são carregadas por contexto, detalhes gênicos por blocos e suporte por região sob demanda. Os downloads maiores usam gzip sem perda de conteúdo. Arquivos antigos substituídos foram retirados da pasta publicada após backup.

- Dados iniciais de RNA: 54,91 → 56,79 MB, acréscimo de aproximadamente 1,88 MB (3,4%).
- Arquivos iniciais do navegador ChIP: aproximadamente 3,18 → 5,14 MB.
- Pacote estático completo: aproximadamente 313 MB, incluindo downloads e detalhes sob demanda; esse total não é transferido na abertura.
- Os números são tamanhos de arquivos, não medições de memória ou latência do navegador.

## Verificação

A validação completa do pacote publicado passou, incluindo checksums, conservação das coordenadas e atribuições nativas, suporte, núcleos, IDR, resumos compactos, genes da referência fora do RNA, links e versões de arquivos. O teste específico de PRJNA312093 conferiu 242.105 regiões e 484.210 registros de suporte. Passaram também os testes de coleções, integração, navegador, trilhas e 11 testes de admissão/conservação. Todos os dados de RNA e sua proveniência foram comparados com a versão anterior e permaneceram iguais.

A revisão visual no localhost conferiu a nova página do estudo, a integração, o detalhe de Smp_243660 em H3K4me3/female_cercaria, cobertura por biblioteca, núcleo compartilhado, IDR separado e restauração das trilhas ao recarregar a URL. O atalho de QC na seleção de genes foi corrigido para respeitar o estudo escolhido; não apareceram erros de execução na aba final. Os hashes finais publicados foram conferidos após essa correção.

Publicação: `results/atlas/rnaseq/`. Backup recuperável: `.codex_tmp/pre_picard_20261021/`. Confirmação visual: `.codex_tmp/picard_atlas_final_20261021.png`. Resultado da validação: `.codex_tmp/picard_validation.json`.
