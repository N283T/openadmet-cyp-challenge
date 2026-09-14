"""Assemble the selected D6 ComplementNB and A4 80/20 rank-blend candidate."""

import json
import os
from pathlib import Path

import polars as pl

from cyp_challenge.tracks import ENDPOINTS, validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


d6 = read("d6_candidate").sort("Molecule_Name")
aliases = (
    read("aliases")
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
a4_scores = (
    read("a4_blends")
    .filter(pl.col("method") == cfg["a4_method"])
    .select("compound_id", "prediction")
)
assert a4_scores.height == 750
assert a4_scores["compound_id"].n_unique() == 750

selected = (
    a4_scores.sort(["prediction", "compound_id"], descending=[True, False])
    .head(cfg["a4_positive_calls"])
    .select("compound_id")
    .with_columns(pl.lit(True).alias("new_a4"))
)
a4_calls = (
    aliases.join(selected, on="compound_id", how="left", validate="1:1")
    .with_columns(pl.col("new_a4").fill_null(False))
    .select("Molecule_Name", "new_a4")
)
candidate = (
    d6.join(a4_calls, on="Molecule_Name", validate="1:1")
    .drop("CYP3A4_is_TDI")
    .rename({"new_a4": "CYP3A4_is_TDI"})
    .select("SMILES", "Molecule_Name", *ENDPOINTS["classification"])
    .sort("Molecule_Name")
)

old_a4 = d6["CYP3A4_is_TDI"].to_numpy().astype(bool)
new_a4 = candidate["CYP3A4_is_TDI"].to_numpy().astype(bool)
summary = {
    "d6_positive_calls": int(candidate["CYP2D6_is_TDI"].sum()),
    "a4_positive_calls": int(candidate["CYP3A4_is_TDI"].sum()),
    "a4_method": cfg["a4_method"],
    "a4_live_only": int((old_a4 & ~new_a4).sum()),
    "a4_candidate_only": int((~old_a4 & new_a4).sum()),
}

candidate.write_parquet(out / "predictions.parquet")
candidate.write_csv(out / "predictions.csv")
summary["validation"] = validate(
    out / "predictions.parquet", "classification", inp / cfg["reference"]
)
summary["uploaded"] = False
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print((out / "summary.json").read_text(), flush=True)
