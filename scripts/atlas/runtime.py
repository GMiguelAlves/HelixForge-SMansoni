"""Small, standard-library-only launcher and static server for the atlas."""
from __future__ import annotations

import argparse
import functools
import hashlib
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from urllib.parse import unquote, urlsplit
from urllib.request import urlopen
import tempfile
import webbrowser

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DIRECTORY = ROOT / "results/atlas/rnaseq"
VERSION = "2026.10.21"
IMAGE = f"helixforge-smansoni-atlas:{VERSION}"
RELEASE = f"https://github.com/GMiguelAlves/HelixForge-SMansoni/releases/download/atlas-v{VERSION}"


def docker_command() -> list[str]:
    docker = shutil.which("docker")
    if docker:
        return [docker]
    wsl = shutil.which("wsl.exe") if os.name == "nt" else None
    if wsl:
        prefix = [wsl, "--exec", "docker"]
        result = subprocess.run(prefix + ["version", "--format", "{{.Client.Version}}"], capture_output=True, timeout=30)
        if result.returncode == 0:
            return prefix
    raise ValueError("Docker was not found. Start Docker Desktop or Docker in the default WSL distribution, or run: python atlas.py init --local")


def load_released_image(docker: list[str]):
    """Download a versioned Docker archive, verify it, then import it locally."""
    filename = f"atlas-docker-{VERSION}.tar.gz"
    with urlopen(RELEASE + "/SHA256SUMS", timeout=60) as response:
        checksums = response.read(16384).decode()
    matches = [line.split()[0] for line in checksums.splitlines()
               if len(line.split()) == 2 and line.split()[1] == filename]
    if len(matches) != 1 or len(matches[0]) != 64:
        raise ValueError("Release checksum is missing or invalid.")
    print("Downloading the frozen atlas Docker image (first launch only)…", flush=True)
    with tempfile.TemporaryDirectory(prefix="hf-atlas-") as temporary:
        archive = Path(temporary) / filename
        digest = hashlib.sha256()
        with urlopen(RELEASE + "/" + filename, timeout=60) as response, archive.open("wb") as output:
            total = 0
            while block := response.read(1024 * 1024):
                total += len(block)
                if total > 2 * 1024**3:
                    raise ValueError("Docker archive exceeds the release size limit.")
                digest.update(block)
                output.write(block)
        if digest.hexdigest() != matches[0]:
            raise ValueError("Docker archive checksum failed; image was not imported.")
        with archive.open("rb") as handle:
            # Streaming also works when Windows Python calls Docker inside WSL.
            subprocess.run(docker + ["load"], stdin=handle, check=True)


def identity(directory: Path) -> str:
    return hashlib.sha256((directory / "manifest.json").read_bytes()).hexdigest()


class AtlasHandler(SimpleHTTPRequestHandler):
    extensions_map = {**SimpleHTTPRequestHandler.extensions_map,
                      ".js": "text/javascript", ".json": "application/json",
                      ".woff2": "font/woff2", ".tsv": "text/tab-separated-values"}

    def __init__(self, *args, directory, **kwargs):
        self.atlas_directory = Path(directory).resolve()
        super().__init__(*args, directory=str(self.atlas_directory), **kwargs)

    def list_directory(self, path):
        self.send_error(404)
        return None

    def send_head(self):
        parsed = urlsplit(self.path)
        if parsed.path == "/healthz":
            data = json.dumps({"application": "helixforge-smansoni-atlas",
                               "manifest_sha256": identity(self.atlas_directory)}).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            import io
            return io.BytesIO(data)
        target = (self.atlas_directory / unquote(parsed.path).lstrip("/")).resolve()
        if (not target.is_relative_to(self.atlas_directory)
                or any(part.startswith(".") for part in Path(unquote(parsed.path)).parts)):
            self.send_error(404)
            return None
        if parsed.path == "/":
            self.path = "/atlas.html" + ("?" + parsed.query if parsed.query else "")
        return super().send_head()

    def end_headers(self):
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()


def health(url: str) -> dict | None:
    try:
        with urlopen(url + "healthz", timeout=2) as response:
            data = json.load(response)
        return data if data.get("application") == "helixforge-smansoni-atlas" else None
    except (OSError, ValueError):
        return None


def open_atlas(url: str, no_browser: bool):
    print(f"Atlas: {url}", flush=True)
    if not no_browser:
        webbrowser.open(url)


def serve(args):
    directory = args.directory.resolve()
    if not (directory / "atlas.html").is_file() or not (directory / "manifest.json").is_file():
        raise ValueError(f"Atlas snapshot is incomplete: {directory}")
    handler = functools.partial(AtlasHandler, directory=directory)
    with ThreadingHTTPServer((args.bind, args.port), handler) as server:
        url = f"http://127.0.0.1:{server.server_port}/"
        open_atlas(url, args.no_browser)
        print("Keep this terminal open. Press Ctrl+C to stop.", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


def init_docker(args):
    docker = docker_command()
    subprocess.run(docker + ["info"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    url = f"http://127.0.0.1:{args.port}/"
    name = f"hf-smansoni-atlas-{args.port}"
    result = subprocess.run(docker + ["container", "inspect", name], capture_output=True, text=True)
    if result.returncode == 0:
        container = json.loads(result.stdout)[0]
        if container.get("Config", {}).get("Labels", {}).get("org.helixforge.application") != "atlas":
            raise ValueError(f"Container {name} is not managed by this launcher.")
        if container.get("Config", {}).get("Image") != args.image:
            raise ValueError(f"Container {name} uses another image. Stop it with: python atlas.py stop --port {args.port}")
        if not container["State"]["Running"]:
            subprocess.run(docker + ["start", name], check=True)
    else:
        present = subprocess.run(docker + ["image", "inspect", args.image], capture_output=True)
        if present.returncode:
            if args.image == IMAGE:
                load_released_image(docker)
            else:
                subprocess.run(docker + ["pull", args.image], check=True)
        subprocess.run(docker + ["run", "--detach", "--name", name,
                        "--label", "org.helixforge.application=atlas",
                        "--read-only", "--cap-drop=ALL", "--security-opt=no-new-privileges",
                        "--publish", f"127.0.0.1:{args.port}:8080", args.image], check=True)
    for _ in range(60):
        current = health(url)
        if current:
            if args.image == IMAGE and (DEFAULT_DIRECTORY / "manifest.json").is_file():
                if current["manifest_sha256"] != identity(DEFAULT_DIRECTORY):
                    raise ValueError("The running container has a different snapshot. Choose a different port or stop the old container.")
            open_atlas(url, args.no_browser)
            return
        time.sleep(0.5)
    raise ValueError(f"Container started, but the atlas did not respond at {url}. Inspect: docker logs {name}")


def stop(args):
    docker = docker_command()
    name = f"hf-smansoni-atlas-{args.port}"
    result = subprocess.run(docker + ["container", "inspect", name], capture_output=True, text=True)
    if result.returncode:
        print("No atlas container at this port.")
        return
    container = json.loads(result.stdout)[0]
    if container.get("Config", {}).get("Labels", {}).get("org.helixforge.application") != "atlas":
        raise ValueError(f"Container {name} is not managed by this launcher.")
    subprocess.run(docker + ["rm", "--force", name], check=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Open the HelixForge-SMansoni atlas.")
    commands = parser.add_subparsers(dest="command", required=True)
    launch = commands.add_parser("init", help="Start Docker and open the atlas; --local uses Python instead.")
    launch.add_argument("--local", action="store_true", help="Serve the included snapshot without Docker.")
    launch.add_argument("--image", default=IMAGE)
    server = commands.add_parser("serve", help="Serve the static snapshot (used inside Docker).")
    server.add_argument("--bind", default="127.0.0.1")
    for command in (launch, server):
        command.add_argument("--port", type=int, default=8765)
        command.add_argument("--directory", type=Path, default=DEFAULT_DIRECTORY)
        command.add_argument("--no-browser", action="store_true")
    shutdown = commands.add_parser("stop", help="Remove the atlas Docker container; keep all data.")
    shutdown.add_argument("--port", type=int, default=8765)
    args = parser.parse_args(argv)
    if not 1 <= args.port <= 65535:
        parser.error("Port must be between 1 and 65535.")
    try:
        if args.command == "stop":
            stop(args)
        elif args.command == "serve" or args.local:
            if args.command == "init":
                args.bind = "127.0.0.1"
                existing = health(f"http://127.0.0.1:{args.port}/")
                if existing:
                    if existing["manifest_sha256"] != identity(args.directory):
                        raise ValueError("Another atlas snapshot is using this port. Choose a different --port.")
                    open_atlas(f"http://127.0.0.1:{args.port}/", args.no_browser)
                    return 0
            serve(args)
        else:
            if args.directory != DEFAULT_DIRECTORY:
                raise ValueError("--directory applies only to --local; Docker uses its frozen snapshot.")
            init_docker(args)
    except (OSError, ValueError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
        print(f"Could not start atlas: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
