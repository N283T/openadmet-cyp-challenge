"""Full-fit D6 one-sided-primary model and current-vector matched candidate."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from regression_d6_primary_upper_bound import fit_model
from scipy.stats import spearmanr


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    scores = pl.read_parquet(inp / "scores.parquet")
    medium = pl.read_parquet(inp / "medium.parquet")
    score_cols = sorted(c for c in scores.columns if c.startswith("score_"))
    medium_cols = sorted(c for c in medium.columns if c.startswith("emb_"))
    cols = score_cols + medium_cols
    if len(cols) != 647:
        raise ValueError("Unexpected frozen features")
    features = scores.select("compound_id", *score_cols).join(
        medium.select("compound_id", *medium_cols), on="compound_id", validate="1:1"
    )
    direct = (
        pl.read_parquet(inp / "measurements.parquet")
        .filter(pl.col("endpoint") == cfg["endpoint"])
        .select("compound_id", "value")
    )
    primary = (
        pl.read_parquet(inp / "regression-labels.parquet")
        .filter(pl.col("head") == cfg["primary_head"])
        .select("compound_id", pl.col("value").alias("primary_value"))
    )
    if direct.height != 1493 or primary.height < 4000:
        raise ValueError("Unexpected D6 source coverage")
    direct_data = direct.join(features, on="compound_id", validate="1:1").sort("compound_id")
    pseudo_data = (
        primary.filter(pl.col("primary_value") >= cfg["primary_weak_cutoff"])
        .join(features, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    if direct_data.join(pseudo_data.select("compound_id"), on="compound_id", how="inner").height:
        raise ValueError("Direct/pseudo compound overlap")
    aliases = (
        pl.read_parquet(inp / "aliases.parquet")
        .filter(pl.col("role") == "test")
        .select("compound_id", "molecule_name")
        .unique()
        .sort("compound_id")
    )
    reference = pl.read_csv(inp / "reference.csv").sort("Molecule_Name")
    current = pl.read_csv(inp / "current-submission.csv").sort("Molecule_Name")
    if (
        reference.height != 750
        or current.height != 750
        or not reference.equals(current.select(reference.columns))
    ):
        raise ValueError("Current submission and official reference differ")
    test = (
        aliases.join(
            current.select("Molecule_Name", "SMILES", cfg["endpoint"]),
            left_on="molecule_name",
            right_on="Molecule_Name",
            validate="1:1",
        )
        .join(features, on="compound_id", validate="1:1")
        .sort("molecule_name")
    )
    if test.height != 750 or test["compound_id"].n_unique() != 750:
        raise ValueError("Invalid blind compound coverage")
    raw = fit_model(
        direct_data.select(cols).to_numpy().astype(np.float32),
        direct_data["value"].to_numpy(),
        pseudo_data.select(cols).to_numpy().astype(np.float32),
        test.select(cols).to_numpy().astype(np.float32),
        cfg["constraint_weight"],
        cfg["direct_upper_bound"],
        cfg["model"],
        cfg["model"]["seed"],
        out / "d6-upper-model.txt",
    )
    old = test[cfg["endpoint"]].to_numpy()
    matched = (raw - raw.mean()) / raw.std() * old.std() + old.mean()
    candidate_d6 = (1 - cfg["blend_weight"]) * old + cfg["blend_weight"] * matched
    replacement = (
        test.select("molecule_name")
        .with_columns(pl.Series(cfg["endpoint"], candidate_d6))
        .rename({"molecule_name": "Molecule_Name"})
    )
    candidate = (
        current.drop(cfg["endpoint"])
        .join(replacement, on="Molecule_Name", validate="1:1")
        .select(current.columns)
        .sort("Molecule_Name")
    )
    for endpoint in current.columns[2:]:
        if endpoint != cfg["endpoint"] and not candidate[endpoint].equals(current[endpoint]):
            raise ValueError(f"Changed unrelated endpoint: {endpoint}")
    if (
        candidate.select("SMILES", "Molecule_Name").equals(
            current.select("SMILES", "Molecule_Name")
        )
        is False
    ):
        raise ValueError("Candidate identity mismatch")
    candidate.write_csv(out / "candidate.csv")
    candidate.write_parquet(out / "candidate.parquet")
    test.select("compound_id", "molecule_name").with_columns(
        pl.lit(cfg["endpoint"]).alias("endpoint"),
        pl.lit("blind").alias("split_id"),
        pl.lit(-1).alias("fold"),
        pl.Series("current_prediction", old),
        pl.Series("raw_prediction", raw),
        pl.Series("matched_prediction", matched),
        pl.Series("prediction", candidate_d6),
        pl.Series("delta", candidate_d6 - old),
    ).write_parquet(out / "d6-test-comparison.parquet")
    summary = {
        "config": cfg,
        "direct_train": direct_data.height,
        "upper_bound_train": pseudo_data.height,
        "test_n": test.height,
        "current_mean": float(old.mean()),
        "current_sd": float(old.std()),
        "raw_mean": float(raw.mean()),
        "raw_sd": float(raw.std()),
        "candidate_mean": float(candidate_d6.mean()),
        "candidate_sd": float(candidate_d6.std()),
        "candidate_vs_current_spearman": float(spearmanr(candidate_d6, old).statistic),
        "mean_absolute_delta": float(np.abs(candidate_d6 - old).mean()),
        "max_absolute_delta": float(np.abs(candidate_d6 - old).max()),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
