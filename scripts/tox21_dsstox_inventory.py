"""Inspect original DSSTox archive fields without discarding provenance."""

import json
import os
import zipfile
from pathlib import Path

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
with zipfile.ZipFile(next(inp.glob("*-dsstox-original.zip"))) as z:
    members = z.namelist()
    (out / "members.json").write_text(json.dumps(members, indent=2))
    print(members, flush=True)
    for n in members:
        if not n.endswith("/"):
            (out / Path(n).name).write_bytes(z.read(n))
            print(n, z.read(n)[:3000], flush=True)
