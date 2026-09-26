"""One-sided weak-primary constraints for direct-pIC50 boosting."""

import json
import os
from pathlib import Path

import lightgbm as lgb
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


def fit_model(
    direct_x: np.ndarray,
    direct_y: np.ndarray,
    pseudo_x: np.ndarray,
    held_x: np.ndarray,
    weight: float,
    upper_bound: float,
    params: dict,
    seed: int,
    model_path: Path | None = None,
) -> np.ndarray:
    center = float(direct_y.mean())
    if weight == 0:
        x = direct_x
        labels = direct_y - center
        is_pseudo = np.zeros(len(labels), dtype=bool)
    else:
        x = np.concatenate([direct_x, pseudo_x])
        labels = np.concatenate(
            [
                direct_y - center,
                np.full(len(pseudo_x), upper_bound - center),
            ]
        )
        is_pseudo = np.concatenate(
            [
                np.zeros(len(direct_y), dtype=bool),
                np.ones(len(pseudo_x), dtype=bool),
            ]
        )
    if len(pseudo_x) < 100:
        raise ValueError("Insufficient weak-primary unlabeled rows")

    def objective(pred: np.ndarray, dataset: lgb.Dataset) -> tuple[np.ndarray, np.ndarray]:
        target = dataset.get_label()
        diff = pred - target
        grad = 2.0 * diff
        hess = np.full(len(diff), 2.0)
        if is_pseudo.any():
            violation = np.maximum(diff[is_pseudo], 0)
            grad[is_pseudo] = 2.0 * weight * violation
            hess[is_pseudo] = np.where(violation > 0, 2.0 * weight, 1e-6)
        return grad, hess

    model = lgb.train(
        {
            "objective": objective,
            "verbosity": -1,
            "num_leaves": params["num_leaves"],
            "learning_rate": params["learning_rate"],
            "min_child_samples": params["min_child_samples"],
            "lambda_l2": params["reg_lambda"],
            "num_threads": params["num_threads"],
            "seed": seed,
            "feature_fraction_seed": seed,
            "bagging_seed": seed,
        },
        lgb.Dataset(x, label=labels, free_raw_data=False),
        num_boost_round=params["n_estimators"],
    )
    if model_path is not None:
        model.save_model(str(model_path))
    return model.predict(held_x) + center


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    names = cfg["endpoints"]
    direct = (
        pl.read_parquet(inp / "measurements.parquet")
        .filter(pl.col("endpoint").is_in([f"{name}_pIC50_direct_inhibition" for name in names]))
        .with_columns(pl.col("endpoint").str.extract(r"^(CYP[^_]+)", 1).alias("name"))
        .select(
            "compound_id",
            "name",
            "endpoint",
            pl.col("value").alias("y_true"),
            "conf_low",
            "conf_high",
        )
    )
    primary = (
        pl.read_parquet(inp / "regression-labels.parquet")
        .filter(pl.col("head").is_in([f"{name}@0.0000495049505" for name in names]))
        .with_columns(pl.col("head").str.extract(r"^(CYP[^@]+)", 1).alias("name"))
        .select("compound_id", "name", pl.col("value").alias("primary_value"))
    )
    scores = pl.read_parquet(inp / "scores.parquet")
    medium = pl.read_parquet(inp / "medium.parquet")
    score_cols = sorted(c for c in scores.columns if c.startswith("score_"))
    medium_cols = sorted(c for c in medium.columns if c.startswith("emb_"))
    cols = score_cols + medium_cols
    features = scores.select("compound_id", *score_cols).join(
        medium.select("compound_id", *medium_cols), on="compound_id", validate="1:1"
    )
    if len(cols) != 647:
        raise ValueError("Unexpected frozen feature dimensions")
    outputs = []
    inventory = []
    for split in cfg["splits"]:
        assignments = pl.read_parquet(inp / f"{split}-assignments.parquet").select(
            "compound_id", "fold", "is_anchor"
        )
        for name in names:
            target_direct = direct.filter(pl.col("name") == name).drop("name")
            target_primary = primary.filter(pl.col("name") == name).drop("name")
            data = (
                features.join(assignments, on="compound_id", validate="1:1")
                .join(target_direct, on="compound_id", how="left", validate="1:1")
                .join(target_primary, on="compound_id", how="left", validate="1:1")
                .sort("compound_id")
            )
            if data.filter(pl.col("y_true").is_not_null()).height < 1200:
                raise ValueError(f"Incomplete direct labels for {name}")
            for fold in range(5):
                fit_mask = (pl.col("fold") != fold) | pl.col("is_anchor")
                held = data.filter(~fit_mask & pl.col("y_true").is_not_null())
                if name == "CYP2D6":
                    direct_fit = data.filter(fit_mask & pl.col("y_true").is_not_null())
                else:
                    direct_fit = data.filter(
                        fit_mask
                        & pl.col("y_true").is_not_null()
                        & (pl.col("primary_value") < cfg["primary_weak_cutoff"])
                    )
                pseudo_fit = data.filter(
                    fit_mask & (pl.col("primary_value") >= cfg["primary_weak_cutoff"])
                )
                if held.height < 100 or direct_fit.height < 500 or pseudo_fit.height < 100:
                    raise ValueError(f"Insufficient {split}/{name}/fold{fold} coverage")
                direct_x = direct_fit.select(cols).to_numpy().astype(np.float32)
                direct_y = direct_fit["y_true"].to_numpy()
                pseudo_x = pseudo_fit.select(cols).to_numpy().astype(np.float32)
                held_x = held.select(cols).to_numpy().astype(np.float32)
                inventory.append(
                    {
                        "split_id": split,
                        "target": name,
                        "fold": fold,
                        "direct_train": direct_fit.height,
                        "pseudo_upper_train": pseudo_fit.height,
                        "held_direct": held.height,
                        "held_weak": held.filter(
                            pl.col("primary_value") >= cfg["primary_weak_cutoff"]
                        ).height,
                    }
                )
                for weight in cfg["constraint_weights"]:
                    prediction = fit_model(
                        direct_x,
                        direct_y,
                        pseudo_x,
                        held_x,
                        weight,
                        cfg["direct_upper_bound"],
                        cfg["model"],
                        cfg["model"]["seed"] + fold,
                    )
                    outputs.append(
                        held.select(
                            "compound_id",
                            "endpoint",
                            "fold",
                            "y_true",
                            "conf_low",
                            "conf_high",
                            "primary_value",
                        ).with_columns(
                            pl.lit(split).alias("split_id"),
                            pl.lit(name).alias("target"),
                            pl.lit(weight).alias("constraint_weight"),
                            pl.Series("prediction", prediction),
                        )
                    )
                    print(split, name, fold, weight, flush=True)
    oof = pl.concat(outputs).sort("split_id", "target", "constraint_weight", "compound_id")
    if oof.select("split_id", "target", "constraint_weight", "compound_id").is_duplicated().any():
        raise ValueError("Duplicate OOF keys")
    rows = []
    for (split, target, weight), frame in oof.group_by("split_id", "target", "constraint_weight"):
        for subset, part in {
            "all": frame,
            "selected_like": frame.filter(pl.col("primary_value") < cfg["primary_weak_cutoff"]),
            "weak_screen": frame.filter(pl.col("primary_value") >= cfg["primary_weak_cutoff"]),
        }.items():
            if part.is_empty():
                continue
            y = part["y_true"].to_numpy()
            pred = part["prediction"].to_numpy()
            lo = part["conf_low"].to_numpy()
            hi = part["conf_high"].to_numpy()
            rows.append(
                {
                    "split_id": split,
                    "target": target,
                    "weight": weight,
                    "subset": subset,
                    "n": part.height,
                    "ST_RAE": st_rae(part) if subset == "all" else None,
                    "MAE": float(np.abs(y - pred).mean()),
                    "soft_error": float(
                        (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).mean()
                    ),
                    "Spearman": float(spearmanr(y, pred).statistic),
                    "mean_prediction": float(pred.mean()),
                }
            )
    oof.write_parquet(out / "oof.parquet")
    pl.DataFrame(rows).sort("split_id", "target", "weight", "subset").write_csv(out / "metrics.csv")
    pl.DataFrame(inventory).sort("split_id", "target", "fold").write_csv(out / "inventory.csv")
    (out / "summary.json").write_text(
        json.dumps({"config": cfg, "oof_rows": oof.height}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
