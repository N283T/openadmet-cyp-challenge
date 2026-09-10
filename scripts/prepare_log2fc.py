"""Prepare compound-keyed auxiliary heads without touching high-fidelity labels."""

import json
import os
from pathlib import Path

import polars as pl

from cyp_challenge.storage import write_json


def main():
    inputs, output = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    raw = pl.read_csv(next(inputs.glob("*-cyp-challenge-single-concentration-TRAIN.csv")))
    aliases = pl.read_parquet(next(inputs.glob("*-aliases.parquet")))
    mapping = (
        aliases.filter(pl.col("role") == "aux")
        .select(
            pl.col("molecule_name").alias("Molecule_Name"),
            pl.col("raw_smiles").alias("SMILES"),
            "compound_id",
        )
        .unique()
    )
    frame = raw.join(mapping, on=["Molecule_Name", "SMILES"], how="left", validate="m:1")
    assert frame["compound_id"].null_count() == 0
    assert frame["log2fc_estimate"].is_finite().all()
    frame = frame.with_columns(
        (pl.col("enzyme") + pl.lit("@") + pl.col("concentration_M").cast(pl.String)).alias("head")
    )
    targets = (
        frame.group_by("compound_id", "head")
        .agg(pl.col("log2fc_estimate").mean().alias("value"), pl.len().alias("source_rows"))
        .sort("compound_id", "head")
    )
    targets.write_parquet(output / "targets.parquet")
    frame.write_parquet(output / "source_rows.parquet")
    heads = (
        targets.group_by("head")
        .agg(
            pl.len().alias("n"),
            pl.col("value").mean(),
            pl.col("source_rows").max().alias("max_replicates"),
        )
        .sort("head")
    )
    write_json(
        output / "summary.json",
        {
            "config": config,
            "source_rows": raw.height,
            "compounds": targets["compound_id"].n_unique(),
            "heads": heads.to_dicts(),
            "replicate_policy": "mean per compound/enzyme/concentration; source rows retained",
        },
    )
    print(heads, flush=True)


if __name__ == "__main__":
    main()
