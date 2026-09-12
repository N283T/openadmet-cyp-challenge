"""Capture pinned split, normalization and dataset implementation support."""

import json
import os
from pathlib import Path

import httpx

out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
rows = []
with httpx.Client(timeout=60, follow_redirects=True) as client:
    for rev in cfg["revisions"]:
        for p in cfg["paths"]:
            url = f"https://raw.githubusercontent.com/{cfg['repository']}/{rev}/{p}"
            r = client.get(url)
            rows.append(dict(revision=rev, path=p, status=r.status_code))
            if r.status_code == 404:
                continue
            r.raise_for_status()
            dest = out / rev / p
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(r.content)
(out / "manifest.json").write_text(json.dumps(rows, indent=2))
