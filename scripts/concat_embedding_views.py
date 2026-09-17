"""Join independently trained embedding views by compound identity."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
joined = None
manifest = []
for view in cfg["views"]:
    frame = pl.read_parquet(inp / view["file"])
    columns = sorted(column for column in frame.columns if column.startswith("emb_"))
    if not columns:
        raise ValueError(f"No embedding columns for {view['name']}")
    if frame["compound_id"].n_unique() != frame.height:
        raise ValueError(f"Duplicate compounds for {view['name']}")
    renamed = {column: f"emb_{view['name']}_{index:04d}" for index, column in enumerate(columns)}
    selected = frame.select("compound_id", *columns).rename(renamed).sort("compound_id")
    if joined is None:
        joined = selected
    else:
        joined = joined.join(selected, on="compound_id", validate="1:1")
    manifest.append(
        {
            "name": view["name"],
            "file": view["file"],
            "input_dimensions": len(columns),
            "output_prefix": f"emb_{view['name']}_",
        }
    )

feature_columns = [column for column in joined.columns if column.startswith("emb_")]
if not np.isfinite(joined.select(feature_columns).to_numpy()).all():
    raise ValueError("Non-finite concatenated embedding")
joined.write_parquet(out / "concatenated-embeddings.parquet")
(out / "manifest.json").write_text(
    json.dumps(
        {
            "views": manifest,
            "compounds": joined.height,
            "dimensions": len(feature_columns),
        },
        indent=2,
    )
    + "\n"
)
print(json.dumps({"compounds": joined.height, "dimensions": len(feature_columns)}))
