# HelixForge-SMansoni Atlas

Explore RNA expression, biological contexts, study designs and chromatin regions on the frozen SM_V10 / WBPS19 reference. No sequencing pipeline or analysis software is needed to consult the atlas.

## Open with Docker

Install Python 3.10+ and start Docker Desktop (Linux containers), Docker Engine on Linux, or Docker in your default WSL distribution. Windows initialization detects Docker inside WSL when a Windows Docker executable is unavailable. Download and extract **atlas-init.zip** from the [atlas release](https://github.com/GMiguelAlves/HelixForge-SMansoni/releases/tag/atlas-v2026.10.21). In its extracted folder, run:

```console
python atlas.py init
```

The first launch downloads the versioned Docker image from the release, checks its SHA-256 and imports it into Docker. It then opens **http://127.0.0.1:8765/**. Subsequent launches reuse the container without downloading the image again. The terminal can be closed after Docker initialization finishes.

No container registry login is required. The prebuilt image targets Linux x86-64; ARM machines can build the Dockerfile locally for their architecture. Stop the container with:

```console
python atlas.py stop
```

This removes only the container managed by the launcher. The image remains available. If port 8765 is occupied, use `python atlas.py init --port 8766` and the same port when stopping. Use `--no-browser` on a headless machine.

## Open locally without Docker

Download and extract **atlas-local.zip** from the same release, or clone this repository. Run:

```console
python atlas.py init --local
```

This opens the same atlas and scientific data at **http://127.0.0.1:8765/**. Keep the terminal open; Ctrl+C stops the server. Local mode requires only Python 3.10+ and works offline after download. Docker works offline after its first successful initialization. Substitute `python3` where that is the Python executable name.

Both modes bind only to this computer by default. The server exposes the static atlas folder, not the repository or personal files. Display preferences and user gene collections are stored in the browser. Neither mode refits models or modifies scientific results.

## Build Docker from the repository

```console
docker compose up --build -d
```

Open **http://127.0.0.1:8765/**. Stop this Compose deployment with `docker compose down`. The image contains only the reviewed atlas, the small Python server and their licenses. It runs as an unprivileged user with a read-only filesystem and no mounted host folders. Alternatively:

```console
docker build -t helixforge-smansoni-atlas:2026.10.21 .
python atlas.py init
```

## Scientific scope

- 12 RNA studies, 268 biological samples, 535 technical runs and 9,914 quantified genes.
- 71 original within-study contrasts. Cross-study expression views are descriptive; no batch-corrected pooled inference is substituted.
- 118 genes in the machinery preset; users can also assemble other collections.
- Two ChIP studies: HP1 (PRJNA1090249), H3K27me3 and H3K4me3 (PRJNA312093); 244,788 regional candidates in total.
- ChIP catalogs, library support, shared cores and optional H3K4me3 IDR evidence remain separate. Original assignments and assessed/unassessed ambiguity are preserved.
- RNA/ChIP integration matches biological contexts across independent experiments. It does not establish sample pairing, causal regulation or absence of a mark in unassigned genes. Differential binding remains **HOLD**.

Study pages and downloadable tables retain QC, source commits and interpretation limits. PRJNA312093 RNA admission remains restricted to 12 approved schistosomulum libraries and three reviewed contrasts. Adding its ChIP catalogs does not admit its excluded adult/cercaria RNA results.

## Size and provenance

The static snapshot is **313,521,726 bytes** (about 314 MB before ZIP/image compression), including its manifest. The largest browser payload is about 56.8 MB. The independent ChIP explorer starts with roughly 5.4 MB of assets and loads large catalogs only as needed. Archive size and initial page-transfer size are different quantities.

For this first release the snapshot is versioned under `results/atlas/rnaseq/`, alongside frozen inputs and build code. This makes a clone self-contained and permits byte-for-byte checks. It increases repository size: future large updates should move immutable snapshots to release assets or object storage instead of accumulating every full snapshot in Git. Temporary build artifacts are excluded.

`manifest.json` records 450 output checksums and 354 frozen input checksums. Packaging refuses missing, modified or unexpected output files. Git attributes preserve checksum-bearing bytes on Windows and Linux. **SHA256SUMS** covers downloadable ZIPs and the Docker archive. The release's image-inspection JSON records the built container and base layers.

For maintenance, see [the build notes](analyses/atlas/README.md) in the repository. Running the atlas needs none of those analysis dependencies. Packaging copies the accepted snapshot; it never refreshes scientific data from a mutable remote branch.
