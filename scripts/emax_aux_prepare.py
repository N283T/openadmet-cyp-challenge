"""Prepare four direct-inhibition Emax auxiliaries from the audited official CSV."""

import hashlib
import json
import os
from pathlib import Path
from urllib.request import urlopen

import polars as pl


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    latest_path, original_path = inp / cfg["latest_file"], inp / cfg["original_file"]
    latest = pl.read_csv(latest_path, infer_schema=False)
    original = pl.read_csv(original_path, infer_schema=False)
    aliases = pl.read_parquet(inp / cfg["aliases_file"])
    assignments = pl.read_parquet(inp / cfg["assignments_file"])
    compounds = pl.read_parquet(inp / cfg["compounds_file"])
    assert latest.columns == original.columns
    original.join(latest, on=latest.columns, how="anti", nulls_equal=True).write_csv(
        out / "removed-source-rows.csv"
    )
    latest.join(original, on=latest.columns, how="anti", nulls_equal=True).write_csv(
        out / "added-source-rows.csv"
    )
    mapping = aliases.select(
        pl.col("molecule_name").alias("Molecule_Name"),
        pl.col("raw_smiles").alias("SMILES"),
        "compound_id",
    ).unique()
    assert mapping.select("Molecule_Name", "SMILES").unique().height == mapping.height
    mapped = latest.join(mapping, on=["Molecule_Name", "SMILES"], how="left", validate="m:1")
    unmapped = mapped.filter(pl.col("compound_id").is_null())
    unmapped.write_csv(out / "unmapped-source-rows.csv")
    assert unmapped.is_empty()
    assert set(mapped["compound_id"]).issubset(set(compounds["compound_id"]))
    assert set(mapped["compound_id"]).issubset(set(assignments["compound_id"]))
    heads = [
        f"{iso}_EmaxVsPosCtrl_direct_inhibition" for iso in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
    ]
    assert all(head in latest.columns for head in heads)
    rows = []
    for head in heads:
        rows.append(
            mapped.select(
                "compound_id",
                "Molecule_Name",
                "SMILES",
                pl.lit(head).alias("head"),
                pl.col(head).alias("raw_value"),
                pl.col(head + "_conf_low").alias("raw_conf_low"),
                pl.col(head + "_conf_high").alias("raw_conf_high"),
            ).with_columns(
                pl.col("raw_value").cast(pl.Float64, strict=False).alias("value"),
                pl.col("raw_conf_low").cast(pl.Float64, strict=False).alias("conf_low"),
                pl.col("raw_conf_high").cast(pl.Float64, strict=False).alias("conf_high"),
                pl.col("raw_value").str.extract(r"^\s*([<>]=?)", 1).alias("censor_relation"),
            )
        )
    long = pl.concat(rows)
    long.write_parquet(out / "source-rows.parquet")
    long.filter(pl.col("censor_relation").is_not_null()).write_parquet(
        out / "censored-rows.parquet"
    )
    usable = long.filter(pl.col("value").is_finite())
    targets = (
        usable.group_by("compound_id", "head")
        .agg(pl.col("value").median())
        .sort("compound_id", "head")
    )
    targets.write_parquet(out / "targets.parquet")
    usable.group_by("compound_id", "head").agg(
        pl.len().alias("source_rows"),
        pl.col("value").n_unique().alias("distinct_values"),
        pl.col("conf_low").min(),
        pl.col("conf_high").max(),
    ).write_parquet(out / "target-provenance.parquet")
    stats = (
        targets.group_by("head")
        .agg(
            pl.len().alias("n"),
            pl.col("value").mean().alias("mean"),
            pl.col("value").std().alias("std"),
            pl.col("value").min().alias("min"),
            pl.col("value").quantile(0.01).alias("q01"),
            pl.col("value").quantile(0.05).alias("q05"),
            pl.col("value").median().alias("median"),
            pl.col("value").quantile(0.95).alias("q95"),
            pl.col("value").quantile(0.99).alias("q99"),
            pl.col("value").max().alias("max"),
            pl.col("value").n_unique().alias("distinct"),
        )
        .sort("head")
    )
    stats.write_csv(out / "distributions.csv")
    targets.join(
        assignments.select("compound_id", "fold", "is_anchor"), on="compound_id", validate="m:1"
    ).group_by("head", "fold", "is_anchor").len().sort("head", "fold", "is_anchor").write_csv(
        out / "fold-label-counts.csv"
    )
    with urlopen(cfg["official_card_url"], timeout=30) as response:
        official = response.read()
    (out / "official-dataset-card.md").write_bytes(official)
    summary = dict(
        config=cfg,
        original_rows=original.height,
        latest_rows=latest.height,
        unique_molecule_names=latest["Molecule_Name"].n_unique(),
        mapped_compounds=mapped["compound_id"].n_unique(),
        duplicate_name_rows=latest.height - latest["Molecule_Name"].n_unique(),
        duplicate_compound_rows=mapped.height - mapped["compound_id"].n_unique(),
        missing_values=long["raw_value"].null_count(),
        nonnumeric_values=long.filter(
            pl.col("raw_value").is_not_null() & pl.col("value").is_null()
        ).height,
        censored_values=long["censor_relation"].is_not_null().sum(),
        nonfinite_numeric=long.filter(
            pl.col("value").is_not_null() & ~pl.col("value").is_finite()
        ).height,
        targets=targets.height,
        target_compounds=targets["compound_id"].n_unique(),
        heads=heads,
        inverted_intervals=long.filter(pl.col("conf_low") > pl.col("conf_high")).height,
        estimates_outside_ci=usable.filter(
            (pl.col("value") < pl.col("conf_low")) | (pl.col("value") > pl.col("conf_high"))
        ).height,
        original_sha256=hashlib.sha256(original_path.read_bytes()).hexdigest(),
        latest_sha256=hashlib.sha256(latest_path.read_bytes()).hexdigest(),
        official_card_sha256=hashlib.sha256(official).hexdigest(),
        transformation="None: preserve original EmaxVsPosCtrl scale and sign",
        duplicate_policy="Median exact finite direct Emax per standardized compound and head",
        censor_policy="Retain censor tokens; do not create point labels from bounds",
        included_conditions=["direct_inhibition"],
        tdi_targets_used=False,
        distributions=stats.to_dicts(),
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
