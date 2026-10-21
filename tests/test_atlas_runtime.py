"""Exercise HTTP delivery, release integrity and repeat Docker initialization."""
import functools
import gzip
import hashlib
import io
import json
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch
from urllib.error import HTTPError
from urllib.request import Request, urlopen
from http.server import ThreadingHTTPServer
from types import SimpleNamespace

from scripts.atlas import runtime
from scripts.atlas.package_release import verify

TEST_TEMP = Path(tempfile.gettempdir())
TEST_TEMP.mkdir(exist_ok=True)

class DeliveryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=TEST_TEMP)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "atlas.html").write_text("<h1>Frozen atlas</h1>", encoding="utf-8")
        (self.root / "manifest.json").write_text('{"status":"PASS"}', encoding="utf-8")
        (self.root / "data").mkdir()
        self.data = gzip.compress(b"gene_id\tvalue\nSmp_1\t1.2\n", mtime=0)
        (self.root / "data/values.tsv.gz").write_bytes(self.data)
        (self.root / ".private").write_text("hidden")
        handler = functools.partial(runtime.AtlasHandler, directory=self.root)
        self.server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.url = f"http://127.0.0.1:{self.server.server_port}/"
        self.addCleanup(self.server.server_close)
        self.addCleanup(self.server.shutdown)

    def test_root_head_health_and_lossless_gzip_download(self):
        with urlopen(self.url) as response:
            self.assertEqual(response.read(), b"<h1>Frozen atlas</h1>")
            self.assertEqual(response.headers["X-Content-Type-Options"], "nosniff")
        with urlopen(Request(self.url, method="HEAD")) as response:
            self.assertEqual(response.read(), b"")
            self.assertEqual(int(response.headers["Content-Length"]), 21)
        with urlopen(self.url + "data/values.tsv.gz") as response:
            self.assertIsNone(response.headers.get("Content-Encoding"))
            self.assertEqual(response.read(), self.data)
        self.assertEqual(runtime.health(self.url)["manifest_sha256"], runtime.identity(self.root))

    def test_no_listing_hidden_files_or_traversal(self):
        for name in ["data/", ".private", "%2eprivate", "%2e%2e/outside", "missing.html"]:
            with self.subTest(name=name), self.assertRaises(HTTPError) as error:
                urlopen(self.url + name)
            self.assertEqual(error.exception.code, 404)

    def test_repeated_local_init_reuses_only_the_same_snapshot(self):
        with patch.object(runtime, "DEFAULT_DIRECTORY", self.root), patch.object(runtime.webbrowser, "open") as browser:
            self.assertEqual(runtime.main(["init", "--local", "--directory", str(self.root), "--port", str(self.server.server_port)]), 0)
            browser.assert_called_once_with(self.url)
        other = self.root / "other"
        other.mkdir()
        (other / "manifest.json").write_text("another snapshot")
        self.assertEqual(runtime.main(["init", "--local", "--directory", str(other), "--port", str(self.server.server_port), "--no-browser"]), 1)


class DockerTests(unittest.TestCase):
    def test_corrupt_release_archive_is_never_loaded(self):
        filename = f"atlas-docker-{runtime.VERSION}.tar.gz"
        checksum = hashlib.sha256(b"original image").hexdigest()
        with patch.object(runtime.tempfile, "tempdir", str(TEST_TEMP)), patch.object(runtime, "urlopen", side_effect=[io.BytesIO(f"{checksum}  {filename}\n".encode()), io.BytesIO(b"corrupt image")]), patch.object(runtime.subprocess, "run") as process:
            with self.assertRaisesRegex(ValueError, "checksum failed"):
                runtime.load_released_image(["docker"])
            process.assert_not_called()

    def test_repeat_docker_init_does_not_create_another_container(self):
        container = {"Config": {"Image": runtime.IMAGE, "Labels": {"org.helixforge.application": "atlas"}}, "State": {"Running": True}}
        def command(args, **kwargs):
            if args[1:] == ["info"]:
                return SimpleNamespace(returncode=0)
            if args[1:3] == ["container", "inspect"]:
                return SimpleNamespace(returncode=0, stdout=json.dumps([container]))
            self.fail(f"Unexpected Docker mutation: {args}")
        with patch.object(runtime.shutil, "which", return_value="docker"), patch.object(runtime.subprocess, "run", side_effect=command), patch.object(runtime, "health", return_value={"manifest_sha256": runtime.identity(runtime.DEFAULT_DIRECTORY)}), patch.object(runtime.webbrowser, "open"):
            self.assertEqual(runtime.main(["init"]), 0)

    def test_stop_refuses_an_unmanaged_container(self):
        with patch.object(runtime.shutil, "which", return_value="docker"), patch.object(runtime.subprocess, "run", return_value=SimpleNamespace(returncode=0, stdout='[{"Config":{"Labels":{}}}]')) as process:
            self.assertEqual(runtime.main(["stop"]), 1)
            self.assertEqual(process.call_count, 1)

    def test_frozen_bundle_verifies_without_analysis_dependencies(self):
        manifest = verify(runtime.DEFAULT_DIRECTORY)
        self.assertEqual(manifest["genes"], 9914)
        self.assertEqual(len(manifest["studies"]), 12)


if __name__ == "__main__":
    unittest.main()
