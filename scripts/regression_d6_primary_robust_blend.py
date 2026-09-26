"""Fixed blending screen for the robust D6 one-sided primary-response model."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr


def st_rae(frame: pl.DataFrame) -> float:
    y = frame["y_true"].to_numpy()
    pred = frame["prediction"].to_numpy()
    lo = frame["conf_low"].to_numpy()
    hi = frame["conf_high"].to_numpy()
    center = y.mean()
    return float(
        (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
        / (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    )


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    upper = (
        pl.read_parquet(inp / "upper-oof.parquet")
        .filter(pl.col("target") == "CYP2D6")
        .select(
            "compound_id",
            "endpoint",
            "fold",
            "split_id",
            "y_true",
            "conf_low",
            "conf_high",
            "robust_delta",
            "constraint_weight",
            pl.col("prediction").alias("upper_prediction"),
        )
    )
    parent = (
        pl.read_parquet(inp / "parent-control.parquet")
        .filter(
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("scenario") == cfg["parent_control_scenario"])
        )
        .select("compound_id", "endpoint", "fold", pl.col("prediction").alias("base_prediction"))
    )
    scaffold = (
        pl.read_parquet(inp / "scaffold-control.parquet")
        .filter(
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("method") == cfg["scaffold_control_method"])
        )
        .select("compound_id", "endpoint", "fold", pl.col("prediction").alias("base_prediction"))
    )
    output = []
    for split, baseline in [("parent", parent), ("scaffold", scaffold)]:
        for delta in cfg["robust_deltas"]:
            for upper_weight in cfg["constraint_weights"]:
                source = upper.filter(
                    (pl.col("split_id") == split)
                    & (pl.col("robust_delta") == delta)
                    & (pl.col("constraint_weight") == upper_weight)
                ).join(baseline, on=["compound_id", "endpoint", "fold"], validate="1:1")
                if source.height != 1457 or source["base_prediction"].null_count():
                    raise ValueError(f"Incomplete {split}/{upper_weight} keys")
                for placement in cfg["placements"]:
                    if placement == "raw":
                        placed = source["upper_prediction"].to_numpy()
                    elif placement == "fold_moment_match":
                        placed = np.empty(source.height)
                        for fold in range(5):
                            mask = source["fold"].to_numpy() == fold
                            new = source["upper_prediction"].to_numpy()[mask]
                            base = source["base_prediction"].to_numpy()[mask]
                            placed[mask] = (new - new.mean()) / new.std() * base.std() + base.mean()
                    else:
                        raise ValueError(placement)
                    for weight in [0.0, *cfg["blend_weights"]]:
                        prediction = (1 - weight) * source[
                            "base_prediction"
                        ].to_numpy() + weight * placed
                        output.append(
                            source.select(
                                "compound_id",
                                "endpoint",
                                "fold",
                                "split_id",
                                "y_true",
                                "conf_low",
                                "conf_high",
                            ).with_columns(
                                pl.lit(delta).alias("robust_delta"),
                                pl.lit(upper_weight).alias("constraint_weight"),
                                pl.lit(placement).alias("placement"),
                                pl.lit(weight).alias("blend_weight"),
                                pl.Series("prediction", prediction),
                            )
                        )
    oof = pl.concat(output).sort(
        "split_id", "robust_delta", "constraint_weight", "placement", "blend_weight", "compound_id"
    )
    if (
        oof.select(
            "split_id",
            "robust_delta",
            "constraint_weight",
            "placement",
            "blend_weight",
            "compound_id",
        )
        .is_duplicated()
        .any()
    ):
        raise ValueError("Duplicate keys")
    metrics = []
    for (split, delta, upper_weight, placement, blend_weight), frame in oof.group_by(
        "split_id", "robust_delta", "constraint_weight", "placement", "blend_weight"
    ):
        y = frame["y_true"].to_numpy()
        pred = frame["prediction"].to_numpy()
        metrics.append(
            {
                "split_id": split,
                "robust_delta": delta,
                "constraint_weight": upper_weight,
                "placement": placement,
                "blend_weight": blend_weight,
                "n": frame.height,
                "ST_RAE": st_rae(frame),
                "MAE": float(np.abs(y - pred).mean()),
                "Spearman": float(spearmanr(y, pred).statistic),
            }
        )
    oof.write_parquet(out / "oof.parquet")
    pl.DataFrame(metrics).sort("split_id", "ST_RAE").write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(
        json.dumps({"config": cfg, "oof_rows": oof.height}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
