# Decisões de curadoria

## PRJEB32839

O ENA expõe 90 BioSamples/180 runs, mas o artigo e o metadata oficial dos autores definem 75 datasets. O pacote faz uma interseção pelos 150 run accessions publicados no arquivo `RNAseq_metadata.tsv` dos autores e exclui 30 runs não pertencentes ao conjunto congelado. Cada amostra tem duas lanes.

## PRJEB14695

O acesso correto é PRJEB14695/ERP016356. Os 138 runs são agrupados em 23 BioSamples, seis runs técnicos por amostra. Os aliases ENA definem tecido, sexo e infecção single-sex versus mixed-sex. Dois registros têm `sample_title` trocado; o pacote usa o `sample_alias`, coerente com o conjunto completo.

## PRJNA597909

Aliases públicos definem maturidade, tratamento e cinco réplicas por grupo. O DE usa uma variável combinada `condition`, pois a API DE 1.0 aceita um fator de interesse e covariáveis aditivas, mas não um termo de interação formal.

## PRJNA602528

Os dez tempos têm uma biblioteca pública cada. O pacote executa até tximport (`run_mode=import`) e bloqueia DESeq2; dez níveis com uma observação cada formariam um desenho saturado.
