# HelixForge-SMansoni Atlas

Explore RNA expression, biological contexts, study designs and chromatin regions on the frozen SM_V10 / WBPS19 reference. No sequencing pipeline or analysis software is needed to consult the atlas.

## Windows step-by-step

You can consult the atlas without Git, R, Conda or any analysis packages. Choose one of these two ways to run it:

| Mode | What to download from the release | Requirements | First launch |
| --- | --- | --- | --- |
| **Local Python** — simplest if you only want to explore | `atlas-local.zip` (about 82 MB) | Python 3.10+ | `py -3 atlas.py init --local` |
| **Docker** — if you already use Docker | `atlas-init.zip` (about 11 KB) | Python 3.10+ and a running Linux Docker engine | `py -3 atlas.py init` |

Both modes open the same atlas. The Docker initializer downloads another **124 MB** on its first launch. The prebuilt image is for Linux x86-64; on a Windows ARM computer, use local Python or build an image for your architecture from the repository.

### 1. Check Python

Open **Windows Terminal**, **PowerShell** or **Command Prompt** from the Start menu and run:

```powershell
py -3 --version
```

You should see `Python 3.x.x`, with version **3.10 or later**. If Python is missing, install it using the [official Python for Windows downloads](https://www.python.org/downloads/windows/), then close and reopen the terminal before checking again. With the current Python Install Manager, you can explicitly install a runtime using `py install 3.12`; this installation command is for that manager, not the older launcher. See [Python's Windows installation guide](https://docs.python.org/3/using/windows.html) if needed.

If `py` is not recognized but `python --version` shows a suitable version, replace `py -3` with `python` in every command below. For example, `py -3 atlas.py init --local` becomes `python atlas.py init --local`.

### 2. Download the correct ZIP

Open the [Atlas 2026.10.21 release](https://github.com/GMiguelAlves/HelixForge-SMansoni/releases/tag/atlas-v2026.10.21), expand **Assets** if necessary, and choose:

- **[atlas-local.zip](https://github.com/GMiguelAlves/HelixForge-SMansoni/releases/download/atlas-v2026.10.21/atlas-local.zip)** for local Python.
- **[atlas-init.zip](https://github.com/GMiguelAlves/HelixForge-SMansoni/releases/download/atlas-v2026.10.21/atlas-init.zip)** for Docker.

Save the ZIP in a folder you can easily find, such as **Downloads**. You only need one of these ZIPs. The automatic **Source code (zip)** download is for the full repository and is not needed for these steps. Docker users do not need to download or manually extract the `.tar.gz` image: `init` handles it.

### 3. Extract the ZIP

In **File Explorer**, find the downloaded ZIP, right-click it and choose **Extract All… / Extrair tudo…**. Choose a destination and click **Extract / Extrair**. This is the standard [Windows ZIP extraction procedure](https://support.microsoft.com/en-us/windows/experience/storage-filemanagement/zip-and-unzip-files).

Open the extracted folders until you reach **HelixForge-SMansoni-Atlas**, which contains **atlas.py**. For example, after extracting the local ZIP into Downloads:

```text
Downloads
└── atlas-local
    └── HelixForge-SMansoni-Atlas   ← open this folder
        ├── atlas.py
        ├── ATLAS.md
        ├── LICENSE
        ├── scripts
        └── results               ← present in the local ZIP
```

The surrounding folder name depends on where you extracted the ZIP. The important check is that **you can see atlas.py in the folder you open**. The small Docker ZIP has no `results` folder; its atlas data are downloaded inside the image.

Run from the extracted folder, not from the ZIP preview. Keep `atlas.py`, `scripts` and the other extracted contents together; do not move only `atlas.py` onto your desktop.

### 4. Open a terminal in the atlas folder

While File Explorer shows the folder containing **atlas.py**:

1. Click the **address bar** at the top of File Explorer, where the folder path appears.
2. Type `powershell` and press **Enter**.
3. A PowerShell window opens in that folder. Run the following check:

```powershell
dir atlas.py
```

It should list `atlas.py`. If it says the file was not found, return to File Explorer and open the inner `HelixForge-SMansoni-Atlas` folder first.

### 5A. Start without Docker

Use this command if you downloaded **atlas-local.zip**:

```powershell
py -3 atlas.py init --local
```

The terminal prints `Atlas: http://127.0.0.1:8765/` and opens your browser. If a browser window does not open automatically, copy that address into your browser's address bar. You can start exploring genes, studies and collections immediately; no data processing is required.

**Keep the terminal open** while using local mode. To stop the atlas, click the terminal and press **Ctrl+C**. To open it again later, repeat steps 4 and 5A. No new download is needed, and this mode works offline after extraction.

### 5B. Start with Docker

Use this path if you downloaded **atlas-init.zip**. You still need Python on Windows to run the initializer.

Start **Docker Desktop** from the Start menu and wait until its engine is running. Use **Linux containers**, not Windows containers; follow the [official Windows Docker setup guide](https://docs.docker.com/desktop/setup/install/windows-install/) if Docker is not installed yet. A Docker Engine already running in your **default WSL distribution** also works. The initializer uses a Windows `docker` executable when one is available, otherwise it tries Docker through WSL.

In the PowerShell window from step 4, run:

```powershell
py -3 atlas.py init
```

The first launch downloads the versioned Docker image (about 124 MB), checks its SHA-256, imports it and starts the atlas. Wait until the terminal prints `Atlas: http://127.0.0.1:8765/`. The browser opens automatically; if it does not, open the printed address manually. No container registry account or login is required.

After initialization finishes, **you may close the terminal**; the Docker container continues running. Keep Docker running while using the atlas. To reopen the atlas later, start Docker if necessary, open PowerShell in the same folder and run `py -3 atlas.py init` again. It reuses the image and container rather than downloading another copy. Subsequent launches work offline once Python, the image and the Docker engine are available.

To stop and remove this atlas container, run:

```powershell
py -3 atlas.py stop
```

The image remains on your computer for future launches.

### Common problems

| What you see | What to do |
| --- | --- |
| `py` is not recognized | Complete step 1. If `python --version` works, use `python` instead of `py -3`. |
| Python opens Microsoft Store, or no runtime is installed | Finish the official Python installation, reopen the terminal and check the version before running the atlas. |
| `can't open file ... atlas.py` | The terminal is in the wrong folder. Repeat step 4 and check with `dir atlas.py`. |
| `Atlas snapshot is incomplete` | Local mode needs the fully extracted `atlas-local.zip`. The small `atlas-init.zip` is for Docker only. |
| Docker is unavailable or its engine is not running | Start Docker Desktop or your existing WSL Docker engine and retry. If using the local ZIP, you can choose `--local` instead. |
| Port 8765 is already in use, or another snapshot is using it | Choose port 8766 using the commands below. |
| The browser cannot reach the atlas | Open the exact URL printed by `init`. In local mode, keep its terminal open; in Docker mode, keep Docker running. |

For another port, keep the same launch mode and add `--port`:

```powershell
# Local Python
py -3 atlas.py init --local --port 8766

# Docker
py -3 atlas.py init --port 8766
```

Then open **http://127.0.0.1:8766/**. Stop local Python with **Ctrl+C**; stop this Docker instance with `py -3 atlas.py stop --port 8766`.

Both modes bind only to this computer by default. The server exposes the static atlas folder, not the repository or personal files. Display preferences and user gene collections are stored in the browser. Neither mode refits models or modifies scientific results.

## Linux and macOS

Download and extract the same ZIP for your chosen mode. Open a terminal in the extracted folder containing `atlas.py`. With Python 3.10+, run `python3 atlas.py init --local` for the local ZIP, or `python3 atlas.py init` for the Docker initializer with a running Linux Docker engine. Docker users stop it with `python3 atlas.py stop`; local users keep the terminal open and stop with **Ctrl+C**. The prebuilt Docker image is Linux x86-64; ARM users can use local Python or build the image from source for their architecture.

Use `--no-browser` on a headless machine and open the printed URL on that same computer. By default, this is a local-only server, not a remotely hosted website.

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
