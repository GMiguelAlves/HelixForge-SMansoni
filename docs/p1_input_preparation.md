# Preparação das entradas RNA-seq P1

## Escopo desta entrega

Esta entrega organiza os sete acessos P1 ainda não processados. Ela congela
metadados, configurações, decisões de contraste e proveniência, mas **não**
executa os estudos e **não** altera o atlas.

| Acesso | Runs | Amostras biológicas | Layout | Modo | Decisão de DE |
|---|---:|---:|---|---|---|
| `PRJNA312093` | 20 | 20 | single | `full` | cinco contrastes fêmea × macho, dentro de estágio |
| `PRJNA494090` | 24 | 24 | paired | `full` | selecionada × controle em P6, P7 e P8 |
| `PRJNA704646` | 24 | 8 | paired | `import` | bloqueado: três bibliotecas técnicas do mesmo pool de RNA |
| `PRJNA1150591` | 15 | 15 | paired | `full` | MCLZ × controle em adultos de 7 semanas; coorte de 4 semanas será ligada a `PRJNA597909` depois |
| `PRJNA1177852` | 8 | 8 | paired | `import` | bloqueado: uma biblioteca por tempo e sem controles paralelos |
| `PRJNA551487` | 10 | 10 | paired | `full` | FPL-64176 × DMSO, condicionado ao gate de QC |
| `PRJNA428470` | 6 | 6 | paired | `full` | 5-azacitidina × controle em fêmeas adultas |

O inventário soma 107 runs, 91 unidades biológicas e 352.540.122.551 bytes de
FASTQ bruto segundo os snapshots ENA congelados em 9 de outubro de 2026.

O artigo de meclonazepam declara que as cinco amostras MCLZ de quatro semanas
devem ser pareadas à coorte de quatro semanas depositada em `PRJNA597909`, já
processada. Esse contraste entre depósitos não é executado nesta preparação:
ele será reconstruído apenas na etapa posterior de integração, preservando a
proveniência de cada estudo.

## Contratos

- `metadata/rnaseq/<ACESSO>/runs.tsv`: URL, MD5, tamanho e associação run → amostra.
- `metadata/rnaseq/<ACESSO>/samples.tsv`: unidade biológica usada na inferência.
- `metadata/rnaseq/<ACESSO>/metadata.csv`: entrada nativa do HelixForge.
- `config/<ACESSO>/de_spec.json`: contrastes permitidos.
- `config/<ACESSO>/de_spec.blocked.json`: bloqueio explícito de inferência.
- `provenance/<ACESSO>/source_snapshots/ena_read_run.tsv`: resposta ENA bruta.
- `provenance/<ACESSO>/input_manifest.json`: hash, seleção e volume bruto.

Os artefatos derivados e snapshots estão congelados. Não há gerador P1 ativo
no repositório; qualquer regeneração requer nova curadoria explícita. Os testes
em `tests/test_p1_preparation.py` verificam contagens, layouts, contrato de
download, replicação dos contrastes e hashes de origem.

## Single-end

`PRJNA312093` declara `library_layout=single`, `fragment_length_mean=200` e
`fragment_length_sd=80`, conforme exige o contrato atual do HelixForge. Esses
valores são um prior provisório de protocolo, não uma medida relatada pelo
artigo. A execução biológica deve revisar a sensibilidade do resultado a esse
prior antes de aceitar o estudo no atlas. O caminho STAR não deve ser usado;
o caminho suportado é FastQC → Trim Galore → merge → Salmon → tximport →
DESeq2.

## Lista de genes congelada

`config/candidates/histone_modifier_gene_set_v2.tsv` preserva a curadoria
histórica e acrescenta as classes importantes ausentes, separadas em writers,
erasers, readers e grupos de contexto. Os cinco IDs históricos antes pendentes
foram resolvidos e auditados em
`provenance/histone_modifier_id_resolution/historical_id_resolution.tsv`.
Os dois PHDs em `REVIEW_ADD` permanecem exploratórios, e os registros em
`EXCLUDED_REVIEW` não entram no núcleo de candidatos.
