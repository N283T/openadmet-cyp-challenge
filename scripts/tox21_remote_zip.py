"""Inspect a public release ZIP using HTTP ranges; retain retrieved segments."""

import io
import json
import os
import zipfile
from pathlib import Path

import httpx

out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


class RemoteFile(io.RawIOBase):
    def __init__(self):
        self.client = httpx.Client(timeout=60, follow_redirects=True)
        resp = self.client.head(cfg["url"])
        self.size = int(resp.headers["content-length"])
        self.pos = 0
        self.n = 0
        (out / "head.json").write_text(json.dumps(dict(resp.headers)))

    def seek(self, offset, whence=0):
        self.pos = (
            offset if whence == 0 else self.pos + offset if whence == 1 else self.size + offset
        )
        return self.pos

    def tell(self):
        return self.pos

    def seekable(self):
        return True

    def read(self, size=-1):
        size = self.size - self.pos if size < 0 else min(size, self.size - self.pos)
        if size == 0:
            return b""
        resp = self.client.get(
            cfg["url"], headers={"Range": f"bytes={self.pos}-{self.pos + size - 1}"}
        )
        resp.raise_for_status()
        assert resp.status_code == 206, resp.status_code
        assert len(resp.content) == size, (len(resp.content), size)
        (out / f"range-{self.n:04d}.bin").write_bytes(resp.content)
        (out / f"range-{self.n:04d}.json").write_text(json.dumps(dict(resp.headers)))
        self.pos += size
        self.n += 1
        return resp.content


with zipfile.ZipFile(RemoteFile()) as archive:
    members = [
        {"name": x.filename, "bytes": x.file_size, "compressed": x.compress_size}
        for x in archive.infolist()
    ]
    (out / "members.json").write_text(json.dumps(members, indent=2))
    print(json.dumps(members), flush=True)
    for name in cfg.get("members", []):
        (out / Path(name).name).write_bytes(archive.read(name))
