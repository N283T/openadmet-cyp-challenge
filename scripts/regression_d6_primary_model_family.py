"""Fit primary-response models with outer-fold compound exclusion."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from catboost import CatBoostRegressor
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor


def make_model(name: str, cfg: dict, fold: int):
    seed = cfg["seed"] + fold
    if name == "lightgbm":
        return LGBMRegressor(
            objective="regression",
            random_state=seed,
            n_jobs=cfg["threads"],
            verbosity=-1,
            **cfg["lightgbm"],
        )
    if name == "catboost":
        return CatBoostRegressor(
            loss_function="RMSE",
            random_seed=seed,
            thread_count=cfg["threads"],
            verbose=False,
            allow_writing_files=False,
            **cfg["catboost"],
        )
    if name == "extra_trees":
        return ExtraTreesRegressor(random_state=seed, n_jobs=cfg["threads"], **cfg["extra_trees"])
    raise ValueError(name)


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    primary = (
        pl.read_parquet(inp / "regression-labels.parquet")
        .filter(pl.col("head") == cfg["primary_head"])
        .select("compound_id", pl.col("value").alias("primary_truth"))
    )
    direct = (
        pl.read_parquet(inp / "measurements.parquet")
        .filter(pl.col("endpoint") == cfg["direct_endpoint"])
        .select("compound_id", pl.col("value").alias("direct_truth"))
    )
    scores = pl.read_parquet(inp / "scores.parquet")
    medium = pl.read_parquet(inp / "medium.parquet")
    score_cols = sorted(c for c in scores.columns if c.startswith("score_"))
    medium_cols = sorted(c for c in medium.columns if c.startswith("emb_"))
    features = scores.select("compound_id", *score_cols).join(
        medium.select("compound_id", *medium_cols), on="compound_id", validate="1:1"
    )
    if (
        primary.height != 4376
        or direct.height != 1493
        or len(score_cols) != 47
        or len(medium_cols) != 600
    ):
        raise ValueError("Unexpected feature or label coverage")
    outputs = []
    for split in cfg["splits"]:
        assignments = pl.read_parquet(inp / f"{split}-assignments.parquet")
        data = (
            assignments.select("compound_id", "fold", "is_anchor")
            .join(features, on="compound_id", validate="1:1")
            .join(primary, on="compound_id", how="left", validate="1:1")
            .join(direct, on="compound_id", how="left", validate="1:1")
            .sort("compound_id")
        )
        if data.filter(pl.col("primary_truth").is_not_null()).height != 4376:
            raise ValueError(f"Incomplete {split} primary labels")
        x = data.select(score_cols + medium_cols).to_numpy().astype(np.float32)
        y = data["primary_truth"].to_numpy()
        fold_ids = data["fold"].to_numpy()
        anchors = data["is_anchor"].to_numpy()
        if not np.isfinite(x).all():
            raise ValueError(f"Nonfinite {split} features")
        for fold in range(5):
            fit = ((fold_ids != fold) | anchors) & np.isfinite(y)
            valid = (fold_ids == fold) & ~anchors & np.isfinite(y)
            if sum(fit) < 3000 or sum(valid) < 500:
                raise ValueError(f"Incomplete {split}/fold{fold} split")
            for name in cfg["models"]:
                model = make_model(name, cfg, fold)
                model.fit(x[fit], y[fit])
                prediction = model.predict(x[valid])
                if not np.isfinite(prediction).all():
                    raise ValueError(f"Nonfinite {split}/{fold}/{name}")
                outputs.append(
                    data.filter(pl.Series(valid))
                    .select("compound_id", "fold", "primary_truth", "direct_truth")
                    .with_columns(
                        pl.lit(split).alias("split_id"),
                        pl.lit(name).alias("method"),
                        pl.Series("prediction", prediction),
                    )
                )
                print(split, fold, name, flush=True)
    oof = pl.concat(outputs).sort("split_id", "method", "compound_id")
    expected = 2 * len(cfg["models"]) * (4376 - 0)
    # Anchors are always training, so OOF coverage is their complement.
    nonanchor_n = oof.filter(
        (pl.col("split_id") == "parent") & (pl.col("method") == cfg["models"][0])
    ).height
    if nonanchor_n >= 4376 or oof.height != 2 * len(cfg["models"]) * nonanchor_n:
        raise ValueError(f"Incomplete OOF coverage: expected upper bound {expected}")
    if oof.select("split_id", "method", "compound_id").is_duplicated().any():
        raise ValueError("Duplicate OOF keys")
    metrics = []
    for (split, method), frame in oof.group_by("split_id", "method"):
        yt = frame["primary_truth"].to_numpy()
        pred = frame["prediction"].to_numpy()
        paired = frame.filter(pl.col("direct_truth").is_not_null())
        metrics.append(
            {
                "split_id": split,
                "method": method,
                "n": frame.height,
                "primary_RMSE": float(np.mean((yt - pred) ** 2) ** 0.5),
                "primary_MAE": float(np.abs(yt - pred).mean()),
                "primary_Spearman": float(spearmanr(yt, pred).statistic),
                "paired_direct_n": paired.height,
                "pred_primary_direct_Spearman": float(
                    spearmanr(paired["direct_truth"], paired["prediction"]).statistic
                ),
            }
        )
    oof.write_parquet(out / "oof.parquet")
    pl.DataFrame(metrics).sort("split_id", "method").write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "oof_rows": oof.height,
                "qualification": "Primary screen only; no direct-model refit",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
