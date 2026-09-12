"""Download author's public datasets and immutable historical producer sources."""

import json
import os
from pathlib import Path

import httpx
import polars as pl

out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
manifest = []
with httpx.Client(timeout=180, follow_redirects=True) as client:
    for name, key in cfg["datasets"].items():
        for extension in ["parquet", "csv"]:
            url = f"https://workbench-public-data.s3.us-west-2.amazonaws.com/{key}.{extension}"
            r = client.get(url)
            if r.status_code in [403, 404]:
                continue
            r.raise_for_status()
            p = out / f"{name}.{extension}"
            p.write_bytes(r.content)
            f = (
                pl.read_parquet(p)
                if extension == "parquet"
                else pl.read_csv(p, infer_schema_length=None)
            )
            manifest.append(
                dict(
                    name=name,
                    url=url,
                    file=p.name,
                    headers=dict(r.headers),
                    rows=f.height,
                    columns=f.columns,
                )
            )
            print(name, f.shape, flush=True)
            break
        else:
            raise ValueError("Public dataset unavailable: " + key)
    for name, url in cfg["sources"].items():
        r = client.get(url)
        r.raise_for_status()
        (out / name).write_bytes(r.content)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
