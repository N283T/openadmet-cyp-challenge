"""Count unique training and blind molecules by official OCNT ID block."""

import json
import os
from pathlib import Path

import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
aliases = pl.read_parquet(inp / cfg["aliases"])
table = (
    aliases.filter(pl.col("role").is_in(["train", "test"]))
    .select(
        "role",
        "compound_id",
        pl.col("molecule_name").str.slice(5, 3).alias("id_block"),
    )
    .unique()
)
if table.filter(pl.col("role") == "test").height != 750:
    raise ValueError("Unexpected blind alias count")
counts = (
    table.group_by("role", "id_block")
    .agg(pl.col("compound_id").n_unique().alias("molecules"))
    .sort("role", "molecules", descending=[False, True])
)
counts.write_csv(out / "id-blocks.csv")
summary = {
    "train_compounds": table.filter(pl.col("role") == "train")["compound_id"].n_unique(),
    "blind_compounds": table.filter(pl.col("role") == "test")["compound_id"].n_unique(),
    "train_in_test_block": table.filter(
        (pl.col("role") == "train") & (pl.col("id_block") == cfg["test_prefix"])
    )["compound_id"].n_unique(),
    "blind_in_test_block": table.filter(
        (pl.col("role") == "test") & (pl.col("id_block") == cfg["test_prefix"])
    )["compound_id"].n_unique(),
    "note": "ID block documents batch provenance, not a chemical label or assay outcome.",
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
