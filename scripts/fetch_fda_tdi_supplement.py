"""Fetch and inventory the FDA CYP3A4 TDI supplementary dataset."""

import hashlib
import json
import os
import zipfile
from pathlib import Path, PurePosixPath

import httpx

out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
archive = out / "fda-cyp3a4-tdi-supplement.zip"
with httpx.Client(follow_redirects=True, timeout=120) as client:
    response = client.get(cfg["url"])
    response.raise_for_status()
archive_md5 = hashlib.md5(response.content, usedforsecurity=False).hexdigest()
if len(response.content) != cfg["expected_bytes"]:
    raise ValueError(f"Unexpected archive size: {len(response.content)} != {cfg['expected_bytes']}")
if archive_md5 != cfg["expected_md5"]:
    raise ValueError(f"Unexpected archive MD5: {archive_md5}")
archive.write_bytes(response.content)

members = []
extract_root = out / "source"
extract_root.mkdir()
with zipfile.ZipFile(archive) as bundle:
    for info in bundle.infolist():
        path = PurePosixPath(info.filename)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError(f"Unsafe archive path: {info.filename}")
        if info.is_dir():
            continue
        content = bundle.read(info)
        destination = extract_root.joinpath(*path.parts)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
        members.append(
            {
                "path": info.filename,
                "bytes": len(content),
                "sha256": hashlib.sha256(content).hexdigest(),
            }
        )

summary = {
    "source": cfg["source"],
    "doi": cfg["doi"],
    "url": cfg["url"],
    "resolved_url": str(response.url),
    "http_status": response.status_code,
    "archive_bytes": len(response.content),
    "archive_md5": archive_md5,
    "archive_sha256": hashlib.sha256(response.content).hexdigest(),
    "members": members,
    "interpretation": (
        "Raw supplementary source only. No endpoint harmonization or training-label "
        "conversion is performed by this run."
    ),
}
(out / "inventory.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
