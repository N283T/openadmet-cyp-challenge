"""Audit D6 direct-label coverage in the existing other-CYP analog proxy."""

import json
import os
from pathlib import Path

import polars as pl


def summarize(name, frame):
    y = frame["value"]
    return {
        "population": name,
        "n": frame.height,
        "mean_pic50": y.mean(),
        "median_pic50": y.median(),
        "below_4": frame.filter(pl.col("value") < 4).height,
        "at_least_5": frame.filter(pl.col("value") >= 5).height,
        "at_least_6": frame.filter(pl.col("value") >= 6).height,
    }


def main():
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    inputs = Path(os.environ["CYP_INPUT_DIR"])
    outputs = Path(os.environ["CYP_OUTPUT_DIR"])
    pairs = pl.read_parquet(inputs / config["analog_pairs"])
    labels = pl.read_parquet(inputs / config["regression_labels"])
    if pairs.height != 750 or pairs.select("anchor_id").n_unique() != 75:
        raise ValueError("Unexpected anchor/analog pair inventory")
    analogs = pairs.select(pl.col("analog_id").alias("compound_id")).unique()
    d6 = labels.filter(pl.col("head") == config["d6_head"])
    if d6["compound_id"].n_unique() != d6.height:
        raise ValueError("Duplicate D6 direct labels")
    selected = d6.join(analogs, on="compound_id", how="inner")
    pl.DataFrame([summarize("all_d6_labels", d6), summarize("pseudo_analogs", selected)]).write_csv(
        outputs / "coverage.csv"
    )
    (outputs / "summary.json").write_text(
        json.dumps(
            {
                "analog_pairs": pairs.height,
                "unique_analogs": analogs.height,
                "d6_direct_labelled_analogs": selected.height,
                "source_analog_run": config["source_analog_run"],
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
