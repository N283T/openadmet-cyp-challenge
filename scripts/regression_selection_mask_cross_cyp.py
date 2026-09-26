"""Mask weak-primary direct labels and test cross-CYP weak-label supervision."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr


def st_rae(y: np.ndarray, lo: np.ndarray, hi: np.ndarray, pred: np.ndarray) -> float:
    center = y.mean()
    numerator = (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
    denominator = (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    return float(numerator / denominator)


def summarize(frame: pl.DataFrame, name: str) -> dict:
    y = frame["y_true"].to_numpy()
    pred = frame["prediction"].to_numpy()
    lo = frame["conf_low"].to_numpy()
    hi = frame["conf_high"].to_numpy()
    soft = np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)
    result = {
        "subset": name,
        "n": frame.height,
        "mean_truth": float(y.mean()),
        "mean_prediction": float(pred.mean()),
        "mean_soft_error": float(soft.mean()),
        "MAE": float(np.abs(y - pred).mean()),
        "Spearman": float(spearmanr(y, pred).statistic) if len(y) > 3 else None,
    }
    result["ST_RAE"] = st_rae(y, lo, hi, pred) if name == "all" else None
    return result


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    direct = (
        pl.read_parquet(inp / "measurements.parquet")
        .filter(
            pl.col("endpoint").is_in(
                [f"{name}_pIC50_direct_inhibition" for name in cfg["endpoints"]]
            )
        )
        .select("compound_id", "endpoint", pl.col("value").alias("y_true"), "conf_low", "conf_high")
    )
    primary = (
        pl.read_parquet(inp / "regression-labels.parquet")
        .filter(pl.col("head").is_in([f"{name}@0.0000495049505" for name in cfg["endpoints"]]))
        .with_columns(pl.col("head").str.extract(r"^(CYP[^@]+)", 1).alias("name"))
        .select("compound_id", "name", pl.col("value").alias("primary_value"))
    )
    labelled = direct.with_columns(
        pl.col("endpoint").str.extract(r"^(CYP[^_]+)", 1).alias("name")
    ).join(primary, on=["compound_id", "name"], how="left", validate="1:1")
    scores = pl.read_parquet(inp / "scores.parquet")
    medium = pl.read_parquet(inp / "medium.parquet")
    score_cols = sorted(c for c in scores.columns if c.startswith("score_"))
    medium_cols = sorted(c for c in medium.columns if c.startswith("emb_"))
    feature_cols = score_cols + medium_cols
    features = scores.select("compound_id", *score_cols).join(
        medium.select("compound_id", *medium_cols), on="compound_id", validate="1:1"
    )
    if len(feature_cols) != 647 or labelled.height != 1412 + 1285 + 2335:
        raise ValueError("Unexpected source coverage")
    output = []
    inventory = []
    cutoff = cfg["primary_weak_cutoff"]
    for split in cfg["splits"]:
        assignments = pl.read_parquet(inp / f"{split}-assignments.parquet").select(
            "compound_id", "fold", "is_anchor"
        )
        data = (
            labelled.join(assignments, on="compound_id", validate="m:1")
            .join(features, on="compound_id", validate="m:1")
            .sort("name", "compound_id")
        )
        for target in cfg["endpoints"]:
            target_data = data.filter(pl.col("name") == target)
            donors = data.filter((pl.col("name") != target) & (pl.col("primary_value") >= cutoff))
            if target_data.is_empty() or donors.is_empty():
                raise ValueError(f"Empty {split}/{target}")
            for fold in range(5):
                target_fit = target_data.filter((pl.col("fold") != fold) | pl.col("is_anchor"))
                target_held = target_data.filter((pl.col("fold") == fold) & ~pl.col("is_anchor"))
                donor_fit = donors.filter((pl.col("fold") != fold) | pl.col("is_anchor"))
                selected_fit = target_fit.filter(pl.col("primary_value") < cutoff)
                if target_held.height < 100 or selected_fit.height < 500:
                    raise ValueError(f"Insufficient {split}/{target}/fold{fold} labels")
                inventory.append(
                    {
                        "split_id": split,
                        "target": target,
                        "fold": fold,
                        "all_target_train": target_fit.height,
                        "selected_target_train": selected_fit.height,
                        "other_cyp_weak_train": donor_fit.height,
                        "held_weak": target_held.filter(pl.col("primary_value") >= cutoff).height,
                    }
                )
                arms = [
                    ("all_labels_oracle", target_fit, None),
                    ("selected_only", selected_fit, None),
                    ("selected_plus_other_weak", selected_fit, donor_fit),
                ]
                for arm, own, other in arms:
                    if other is None:
                        fit_x = own.select(feature_cols).to_numpy().astype(np.float32)
                        fit_y = own["y_true"].to_numpy()
                        weights = np.ones(own.height)
                    else:
                        fit_x = np.concatenate(
                            [
                                own.select(feature_cols).to_numpy().astype(np.float32),
                                other.select(feature_cols).to_numpy().astype(np.float32),
                            ]
                        )
                        fit_y = np.concatenate(
                            [own["y_true"].to_numpy(), other["y_true"].to_numpy()]
                        )
                        weights = np.concatenate(
                            [np.ones(own.height), np.full(other.height, cfg["donor_weight"])]
                        )
                    params = cfg["model"]
                    model = LGBMRegressor(
                        objective="regression_l1",
                        n_estimators=params["n_estimators"],
                        num_leaves=params["num_leaves"],
                        learning_rate=params["learning_rate"],
                        min_child_samples=params["min_child_samples"],
                        reg_lambda=params["reg_lambda"],
                        n_jobs=params["n_jobs"],
                        random_state=params["seed"] + fold,
                        verbosity=-1,
                    )
                    model.fit(fit_x, fit_y, sample_weight=weights)
                    pred = model.predict(
                        target_held.select(feature_cols).to_numpy().astype(np.float32)
                    )
                    output.append(
                        target_held.select(
                            "compound_id",
                            "endpoint",
                            "fold",
                            "y_true",
                            "conf_low",
                            "conf_high",
                            "primary_value",
                        ).with_columns(
                            pl.lit(split).alias("split_id"),
                            pl.lit(target).alias("target"),
                            pl.lit(arm).alias("arm"),
                            pl.Series("prediction", pred),
                        )
                    )
                    print(split, target, fold, arm, flush=True)
    oof = pl.concat(output).sort("split_id", "target", "arm", "compound_id")
    if oof.select("split_id", "target", "arm", "compound_id").is_duplicated().any():
        raise ValueError("Duplicate OOF keys")
    rows = []
    for (split, target, arm), frame in oof.group_by("split_id", "target", "arm"):
        subsets = {
            "all": frame,
            "selected_like": frame.filter(pl.col("primary_value") < cutoff),
            "weak_screen": frame.filter(pl.col("primary_value") >= cutoff),
            "primary_missing": frame.filter(pl.col("primary_value").is_null()),
        }
        for subset, part in subsets.items():
            if part.is_empty():
                continue
            rows.append(
                {"split_id": split, "target": target, "arm": arm, **summarize(part, subset)}
            )
    oof.write_parquet(out / "oof.parquet")
    pl.DataFrame(rows).sort("split_id", "target", "subset", "arm").write_csv(out / "metrics.csv")
    pl.DataFrame(inventory).sort("split_id", "target", "fold").write_csv(out / "inventory.csv")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "oof_rows": oof.height,
                "scope": (
                    "Synthetic target-label masking; target weak labels never train selected arms"
                ),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
