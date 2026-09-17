"""Create a leak-safe external observation subset for focused encoder training."""

import json
import os
from pathlib import Path

import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    exact = inp / name
    paths = [exact] if exact.exists() else [
        path for path in inp.iterdir() if path.name.partition("-")[2] == name
    ]
    if len(paths) != 1:
        raise ValueError(f"Expected one {name}; got {paths}")
    return pl.read_parquet(paths[0])


compounds = read("external-compounds.parquet")
observations = read("observations.parquet")
selected = observations.filter(pl.col("isoform").is_in(cfg["isoforms"]))
if cfg.get("sources"):
    selected = selected.filter(pl.col("source").is_in(cfg["sources"]))
assert selected.height and selected["observation_id"].n_unique() == selected.height
assert not selected["challenge_overlap"].any()
selected_compounds = compounds.join(
    selected.select("compound_id").unique(), on="compound_id", how="semi"
)
assert selected.join(selected_compounds, on="compound_id", how="anti").height == 0
selected_compounds.write_parquet(out / "external-compounds.parquet")
selected.write_parquet(out / "observations.parquet")
inventory = selected.group_by("source", "isoform", "measurement_type").agg(
    pl.len().alias("observations"),
    pl.col("compound_id").n_unique().alias("compounds"),
    pl.col("assay_id").n_unique().alias("assays"),
    pl.col("value").is_not_null().sum().alias("continuous"),
    pl.col("binary_label").is_not_null().sum().alias("binary"),
).sort("source", "isoform", "measurement_type")
inventory.write_csv(out / "inventory.csv")
summary = {
    "config": cfg,
    "base_compounds": compounds.height,
    "selected_compounds": selected_compounds.height,
    "selected_observations": selected.height,
    "selected_assays": selected["assay_id"].n_unique(),
    "inventory": inventory.to_dicts(),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
