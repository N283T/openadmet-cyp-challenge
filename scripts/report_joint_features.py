"""Join measured-substrate predictions, ligand geometry and Jazzy without labels."""

import json
import os
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
joined, columns = None, {}
for index, prefix in enumerate(["geometry", "substrate", "jazzy"]):
    t = pl.read_parquet(next(inp.glob(f"{index:03d}-*")))
    keep = [c for c in t.columns if c not in ["compound_id", "jazzy_success_count"]]
    t = t.select("compound_id", *keep).rename({c: prefix + "__" + c for c in keep})
    assert t.height == t["compound_id"].n_unique() == 6896
    if joined is not None:
        assert set(t["compound_id"]) == set(joined["compound_id"])
    joined = t if joined is None else joined.join(t, on="compound_id", validate="1:1")
    columns[prefix] = [prefix + "__" + c for c in keep]
joined.sort("compound_id").write_parquet(out / "joint-features.parquet")
(out / "feature-columns.json").write_text(json.dumps(columns, indent=2))
