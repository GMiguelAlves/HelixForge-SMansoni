# Atlas build and admission

The atlas is a static consultation layer over accepted RNA and ChIP results. See [ATLAS.md](../../ATLAS.md) for end-user installation and Docker/local commands.

`atlas_config.json` declares 12 RNA studies and two frozen ChIP snapshots. HP1 is pinned to `c9e77fee805d73e0527a02f83217c7a2763d270c`; PRJNA312093 H3K27me3/H3K4me3 is pinned to `e583c84ef5cd37edd2f4df8450f17cc42c39c279`. Compact source publications are retained under `local_inputs/hf_chipseq/`, with manifests and checksums. The historical atomic HP1 snapshot is retained solely as a regression fixture, outside the active catalog.

Curated RNA inputs are under `config/<study>/`, `metadata/<study>/` and `results/<study>/`. They preserve accepted values, frozen DE contrasts and QC summaries. Literature supplements remain separate from measured bulk RNA expression. Current-reference transcript models and chromosome lengths are under `config/reference/`.

See [PRJNA312093 incorporation](PRJNA312093_incorporation_2026-10-21.md) for its scientific admission. Source ambiguity and QC limitations remain visible; no atlas command refits a model or establishes causal regulation or differential binding.

## Build and validate

Use Python 3.12 and Node.js 22+ for maintenance. From the repository root:

```console
python -m pip install -r requirements-atlas-build.txt
python scripts/atlas/build_rnaseq_atlas.py --root . --config analyses/atlas/atlas_config.json --output .atlas/rebuild
python scripts/atlas/validate_rnaseq_atlas.py . .atlas/rebuild
```

Use a new output directory: the builder replaces the directory passed through `--output`. Validate and review it before replacing `results/atlas/rnaseq/`. Maintainable UI sources are in `rnaseq_atlas_template.html`, `collections/`, `chipseq_portal/` and `design/`. `design/README.md` documents visual decisions and themes.

For presentation-only changes, `scripts/atlas/refresh_atlas_design.py` preserves and checks scientific asset hashes; consult its `--help`. Keep the manifest synchronized with rendered assets.

```console
python scripts/atlas/validate_rnaseq_atlas.py
python scripts/atlas/package_release.py --verify-only
python -m unittest discover -s tests -p 'test_atlas*.py'
node tests/test_atlas_browser_logic.cjs
node tests/test_picard_atlas.cjs
```

Additional JavaScript tests cover collection calculations, biological joins, track grouping, the BigWig parser and themes. The release workflow runs these before building and testing the Docker image.

## Publish another release

Update the version in the launcher, Dockerfile, Compose file and user documentation. Keep science, UI and packaging changes in separate descriptive commits. Verify the snapshot, push commits and push an `atlas-v<version>` tag. The workflow creates ZIPs, builds/tests a Linux x86-64 image and publishes its archive with checksums using the repository's scoped GitHub Actions token. No personal token or container registry secret is required.

Do not replace a published snapshot under an existing tag. Git identifies exact code and input/output bytes; the image-inspection JSON also records the built image and its base layers.

