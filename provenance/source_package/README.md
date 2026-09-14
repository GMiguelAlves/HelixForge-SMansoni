# Frozen source package

`frozen/` is the byte-preserved curated input package created on 11 September
2026 and validated against HelixForge commit
`8a1eac201d34861a9f8d3cd595b472a5ca76cc97`.

Its [`MANIFEST.sha256`](frozen/MANIFEST.sha256) records 51 original relative
paths and content hashes. The active project files were reorganized into
`config/`, `metadata/`, and `scripts/`; the original package remains available
here as immutable provenance rather than being duplicated beside `frozen/`.

Validate it from this repository root with:

```bash
python provenance/source_package/frozen/scripts/validate_package.py \
  provenance/source_package/frozen
```
