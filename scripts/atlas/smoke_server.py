#!/usr/bin/env python3
"""Check the running server against exact frozen RNA/ChIP snapshot bytes."""
import hashlib
import json
from pathlib import Path
import sys
import time
from urllib.error import URLError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[2]
BUNDLE = ROOT / "results/atlas/rnaseq"


def check(url):
    for attempt in range(30):
        try:
            with urlopen(url + "healthz", timeout=3) as response:
                health = json.load(response)
            break
        except (OSError, ValueError):
            if attempt == 29:
                raise
            time.sleep(1)
    manifest_bytes = (BUNDLE / "manifest.json").read_bytes()
    assert health["manifest_sha256"] == hashlib.sha256(manifest_bytes).hexdigest()
    records = {record["path"]: record for record in json.loads(manifest_bytes)["outputs"]}
    for name in ["atlas.html", "chipseq.html", "study_PRJEB32839.html", "chipseq_PRJNA312093.html",
                 "atlas_theme.js", "data/atlas_payload.js", "data/chipseq/explorer.json"]:
        with urlopen(url + name, timeout=60) as response:
            assert hashlib.sha256(response.read()).hexdigest() == records[name]["sha256"], name
    with urlopen(url, timeout=10) as response:
        assert hashlib.sha256(response.read()).hexdigest() == records["atlas.html"]["sha256"]
    gz = next(name for name in records if name.endswith(".tsv.gz"))
    with urlopen(url + gz, timeout=60) as response:
        assert response.headers.get("Content-Encoding") is None
        assert hashlib.sha256(response.read()).hexdigest() == records[gz]["sha256"]
    with urlopen(Request(url + "chipseq.html", method="HEAD"), timeout=10) as response:
        assert response.headers.get_content_type() == "text/html"
    print("PASS: exact RNA, ChIP, study, theme and gzip assets served from the frozen snapshot.")


if __name__ == "__main__":
    check(sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8765/")
