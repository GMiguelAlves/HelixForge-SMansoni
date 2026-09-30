# RNA-seq atlas

The RNA-seq atlas is a descriptive, navigable integration of the seven accepted
project studies. It does not pool studies for differential-expression
inference and does not perform batch correction.

## Contract

- biological samples are the unit of observation;
- all studies must use the frozen SM_V10/WBPS19 gene universe;
- TPM is used for cross-study exploratory displays as `log2(TPM + 1)`;
- differential-expression statistics retain their original per-study DESeq2
  models and contrast orientation;
- `not tested`, `not significant`, and `gene absent` are distinct states;
- global PCA is descriptive and must not be interpreted as a batch test;
- canonical `Smp_*` identifiers are accompanied by WBPS19 functional
  descriptions, biotypes, source accessions and previous stable IDs when
  available; descriptions are not presented as official gene symbols;
- candidate IDs that are absent from the current matrix are resolved only when
  WBPS19 provides exactly one previous-stable-ID mapping; ambiguous mappings
  remain explicit and are never selected automatically;
- server paths, user identifiers, and private operational logs are forbidden
  from the public bundle.

The generated bundle is written to `results/atlas/rnaseq/` and contains a
navigable HTML report, exportable SVG figures, machine-readable tables, a
compact browser payload, and a checksum manifest.

## Build

Create the environment once:

```bash
conda env create -f analyses/atlas/environment.yml
conda activate helixforge-smansoni-atlas
```

Then run:

```bash
python scripts/atlas/build_rnaseq_atlas.py \
  --root . \
  --config analyses/atlas/atlas_config.json \
  --output results/atlas/rnaseq
```

The build is deterministic for unchanged inputs. NumPy is used only for
descriptive PCA and correlation calculations. HTML and SVG generation use the
Python standard library.

## Deferred work

Batch-effect assessment is deliberately outside this version. The atlas keeps
the relevant metadata fields so a later reviewed analysis can quantify
technical and biological associations without rebuilding the inventory.
Pathway enrichment is also deferred until an annotation contract is frozen.
