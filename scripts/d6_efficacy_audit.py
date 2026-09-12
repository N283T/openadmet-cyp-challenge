"""Snapshot community efficacy recipe and inspect preserved AID1851 fields."""

import json
import os
from pathlib import Path

import httpx
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
with httpx.Client(timeout=90, follow_redirects=True) as client:
    for name, url in cfg["urls"].items():
        r = client.get(url)
        r.raise_for_status()
        (out / name).write_bytes(r.content)
for i in range(3):
    p = next(inp.glob(f"{i:03d}-*"))
    f = pl.read_parquet(p)
    (out / f"schema-{i}.json").write_text(
        json.dumps(
            dict(
                file=p.name,
                shape=f.shape,
                schema={k: str(v) for k, v in f.schema.items()},
                sample=f.head(2).to_dicts(),
            ),
            indent=2,
            default=str,
        )
    )
    print(i, f.shape, f.columns, flush=True)
