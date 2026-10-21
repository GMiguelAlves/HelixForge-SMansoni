Explore RNA expression, chromatin regions and user-defined gene collections on SM_V10 / WBPS19. Includes the unified design and light/dark/system themes.

**Docker:** extract `atlas-init.zip`, start Docker and run `python atlas.py init`. The first launch downloads and verifies the image, then opens the atlas. Python 3.10+ and Docker with Linux containers are required; the prebuilt image targets Linux x86-64.

**Local:** extract `atlas-local.zip`, then run `python atlas.py init --local`. Keep the terminal open. No analysis dependencies or Docker are required.

Both open **http://127.0.0.1:8765/**. Use `--port 8766` if necessary. Stop the Docker container with `python atlas.py stop`.

Includes 12 RNA studies, 9,914 quantified genes and 244,788 ChIP regional candidates across HP1, H3K27me3 and H3K4me3. Integration is descriptive across independent experiments; differential binding remains HOLD. The snapshot is about 314 MB before ZIP/image compression. QC, provenance, restricted RNA admissions and gene-assignment ambiguity are retained.

`SHA256SUMS` records artifact checksums. `docker-image.json` identifies the built image and its layers. The image can also be imported manually with `docker load --input atlas-docker-2026.10.21.tar.gz`.

See the [current Windows step-by-step guide](https://github.com/GMiguelAlves/HelixForge-SMansoni/blob/master/ATLAS.md#windows-step-by-step) for downloading, extracting the ZIP, opening PowerShell and starting the atlas with Python or Docker/WSL. The [versioned guide](https://github.com/GMiguelAlves/HelixForge-SMansoni/blob/atlas-v2026.10.21/ATLAS.md) records the original release instructions, scientific scope and size notes.
