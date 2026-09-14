"""Identify the submitted classification row across pre/post LB snapshots."""

import json
import os
from pathlib import Path

import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

specs = {
    "MA": ("MA-MCC", cfg["previous_scores"]["MA"]),
    "CYP2D6_is_TDI": ("MCC", cfg["previous_scores"]["CYP2D6_is_TDI"]),
    "CYP3A4_is_TDI": ("MCC", cfg["previous_scores"]["CYP3A4_is_TDI"]),
}
result = {}
username = None
for endpoint, (metric, previous_score) in specs.items():
    pre = pl.read_csv(inp / cfg["pre"][endpoint])
    post = pl.read_csv(inp / cfg["post"][endpoint])
    matched = pre.filter((pl.col(metric) - previous_score).abs() < 1e-10)
    assert matched.height == 1, (endpoint, matched.height)
    row = matched.row(0, named=True)
    if username is None:
        username = row["Username"]
    assert row["Username"] == username
    updated = post.filter(pl.col("Username") == username)
    assert updated.height == 1
    new_row = updated.row(0, named=True)
    result[endpoint] = {
        "pre_rank": row["Rank"],
        "pre_score": row[metric],
        "post_rank": new_row["Rank"],
        "post_score": new_row[metric],
        "post_submitted": new_row["Submitted"],
        "changed": new_row[metric] != row[metric],
    }
(out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result), flush=True)

