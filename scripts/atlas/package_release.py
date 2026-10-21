#!/usr/bin/env python3
"""Verify and package the reviewed static snapshot, without rerunning analyses."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[2]


def verify(bundle):
    bundle = Path(bundle).resolve()
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    expected = {"manifest.json"}
    if manifest.get("status") != "PASS":
        raise ValueError("Only a validated atlas can be packaged.")
    for record in manifest["outputs"]:
        name = record["path"]
        path = (bundle / name).resolve()
        if not path.is_relative_to(bundle) or name in expected:
            raise ValueError(f"Unsafe or duplicate manifest path: {name}")
        expected.add(name)
        if not path.is_file() or path.stat().st_size != record["bytes"]:
            raise ValueError(f"Frozen snapshot checksum mismatch: {name}")
        with path.open("rb") as handle:
            if hashlib.file_digest(handle, "sha256").hexdigest() != record["sha256"]:
                raise ValueError(f"Frozen snapshot checksum mismatch: {name}")
    actual = {path.relative_to(bundle).as_posix() for path in bundle.rglob("*") if path.is_file()}
    if actual != expected:
        raise ValueError(f"Unexpected or missing snapshot files: {sorted(actual ^ expected)}")
    return manifest


def package(output, root=ROOT):
    bundle = root / "results/atlas/rnaseq"
    manifest = verify(bundle)
    output.mkdir(parents=True, exist_ok=True)
    runtime = ["atlas.py", "scripts/atlas/runtime.py", "LICENSE", "ATLAS.md"]
    for filename, include_data in [("atlas-init.zip", False), ("atlas-local.zip", True)]:
        paths = [root / name for name in runtime]
        if include_data:
            paths += [path for path in bundle.rglob("*") if path.is_file()]
        with zipfile.ZipFile(output / filename, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
            for path in sorted(paths):
                info = zipfile.ZipInfo("HelixForge-SMansoni-Atlas/" + path.relative_to(root).as_posix(), (2026, 10, 21, 0, 0, 0))
                info.external_attr = 0o100644 << 16
                info.compress_type = zipfile.ZIP_DEFLATED
                archive.writestr(info, path.read_bytes())
        print(f"{filename}: {(output / filename).stat().st_size:,} bytes")
    print(f"Verified {len(manifest['outputs'])} frozen outputs; no scientific files changed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / ".atlas/release")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    if args.verify_only:
        verify(ROOT / "results/atlas/rnaseq")
        print("Frozen atlas checksums: PASS")
    else:
        package(args.output.resolve())
