"""Test whether other-CYP weak-direct labels help CYP2D6 held-out labels."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    names = [cfg["target"], *cfg["donors"]]
    direct = (
        pl.read_parquet(inp / "measurements.parquet")
        .filter(pl.col("endpoint").is_in([f"{name}_pIC50_direct_inhibition" for name in names]))
        .select("compound_id", "endpoint", pl.col("value").alias("y_true"), "conf_low", "conf_high")
        .with_columns(pl.col("endpoint").str.extract(r"^(CYP[^_]+)", 1).alias("name"))
    )
    primary = (
        pl.read_parquet(inp / "regression-labels.parquet")
        .filter(pl.col("head").is_in([f"{name}@0.0000495049505" for name in names]))
        .with_columns(pl.col("head").str.extract(r"^(CYP[^@]+)", 1).alias("name"))
        .select("compound_id", "name", pl.col("value").alias("primary_value"))
    )
    labelled = direct.join(primary, on=["compound_id", "name"], how="left", validate="1:1")
    scores = pl.read_parquet(inp / "scores.parquet")
    medium = pl.read_parquet(inp / "medium.parquet")
    score_cols = sorted(c for c in scores.columns if c.startswith("score_"))
    medium_cols = sorted(c for c in medium.columns if c.startswith("emb_"))
    cols = score_cols + medium_cols
    features = scores.select("compound_id", *score_cols).join(
        medium.select("compound_id", *medium_cols), on="compound_id", validate="1:1"
    )
    if len(cols) != 647 or labelled.filter(pl.col("name") == cfg["target"]).height != 1493:
        raise ValueError("Unexpected source coverage")
    cutoff = cfg["primary_weak_cutoff"]
    outputs = []
    inventory = []
    for split in cfg["splits"]:
        assignments = pl.read_parquet(inp / f"{split}-assignments.parquet").select(
            "compound_id", "fold", "is_anchor"
        )
        data = labelled.join(assignments, on="compound_id", validate="m:1").join(
            features, on="compound_id", validate="m:1"
        )
        target_data = data.filter(pl.col("name") == cfg["target"])
        donors = data.filter(
            pl.col("name").is_in(cfg["donors"]) & (pl.col("primary_value") >= cutoff)
        )
        for fold in range(5):
            own = target_data.filter((pl.col("fold") != fold) | pl.col("is_anchor"))
            held = target_data.filter((pl.col("fold") == fold) & ~pl.col("is_anchor"))
            other = donors.filter((pl.col("fold") != fold) | pl.col("is_anchor"))
            if held.height < 100 or other.height < 500:
                raise ValueError(f"Insufficient {split}/fold{fold} data")
            inventory.append(
                {
                    "split_id": split,
                    "fold": fold,
                    "d6_train": own.height,
                    "d6_held": held.height,
                    "other_weak_train": other.height,
                }
            )
            own_x = own.select(cols).to_numpy().astype(np.float32)
            other_x = other.select(cols).to_numpy().astype(np.float32)
            held_x = held.select(cols).to_numpy().astype(np.float32)
            if cfg.get("donor_task_feature", False):
                zeros_own = np.zeros((own.height, len(cfg["donors"])), dtype=np.float32)
                zeros_held = np.zeros((held.height, len(cfg["donors"])), dtype=np.float32)
                donor_tags = np.column_stack(
                    [
                        (other["name"].to_numpy() == donor).astype(np.float32)
                        for donor in cfg["donors"]
                    ]
                )
                own_x = np.concatenate([own_x, zeros_own], axis=1)
                other_x = np.concatenate([other_x, donor_tags], axis=1)
                held_x = np.concatenate([held_x, zeros_held], axis=1)
            for weight in cfg["donor_weights"]:
                if weight == 0:
                    fit_x = own_x
                    fit_y = own["y_true"].to_numpy()
                    weights = None
                else:
                    fit_x = np.concatenate([own_x, other_x])
                    fit_y = np.concatenate([own["y_true"].to_numpy(), other["y_true"].to_numpy()])
                    weights = np.concatenate([np.ones(own.height), np.full(other.height, weight)])
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
                pred = model.predict(held_x)
                outputs.append(
                    held.select(
                        "compound_id", "endpoint", "fold", "y_true", "conf_low", "conf_high"
                    ).with_columns(
                        pl.lit(split).alias("split_id"),
                        pl.lit(weight).alias("donor_weight"),
                        pl.Series("prediction", pred),
                    )
                )
                print(split, fold, weight, flush=True)
    oof = pl.concat(outputs).sort("split_id", "donor_weight", "compound_id")
    if oof.select("split_id", "donor_weight", "compound_id").is_duplicated().any():
        raise ValueError("Duplicate OOF keys")
    metrics = []
    for (split, weight), frame in oof.group_by("split_id", "donor_weight"):
        y = frame["y_true"].to_numpy()
        pred = frame["prediction"].to_numpy()
        lo = frame["conf_low"].to_numpy()
        hi = frame["conf_high"].to_numpy()
        center = y.mean()
        soft = np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)
        baseline = np.maximum(lo - center, 0) + np.maximum(center - hi, 0)
        metrics.append(
            {
                "split_id": split,
                "donor_weight": weight,
                "n": frame.height,
                "ST_RAE": float(soft.sum() / baseline.sum()),
                "mean_soft_error": float(soft.mean()),
                "MAE": float(np.abs(y - pred).mean()),
                "Spearman": float(spearmanr(y, pred).statistic),
                "mean_prediction": float(pred.mean()),
            }
        )
    oof.write_parquet(out / "oof.parquet")
    pl.DataFrame(metrics).sort("split_id", "donor_weight").write_csv(out / "metrics.csv")
    pl.DataFrame(inventory).sort("split_id", "fold").write_csv(out / "inventory.csv")
    (out / "summary.json").write_text(
        json.dumps({"config": cfg, "oof_rows": oof.height}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
