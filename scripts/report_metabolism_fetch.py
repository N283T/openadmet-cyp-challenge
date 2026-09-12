"""Snapshot official TDC substrate data and source/license metadata."""

import json
import os
from pathlib import Path

import httpx

out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
with httpx.Client(follow_redirects=True, timeout=120) as client:
    for name, url in cfg["sources"].items():
        r = client.get(url)
        r.raise_for_status()
        (out / name).write_bytes(r.content)
        print(name, len(r.content), r.text[:150], flush=True)
