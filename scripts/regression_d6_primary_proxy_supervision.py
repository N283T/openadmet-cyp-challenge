"""Use measured primary response as weak supervision for unlabelled direct potency."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LinearRegression


def st_rae(y: np.ndarray, lo: np.ndarray, hi: np.ndarray, pred: np.ndarray) -> float:
    center = y.mean()
    numerator = (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
    denominator = (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    return float(numerator / denominator)


def map_primary(
    name: str, x_train: np.ndarray, y_train: np.ndarray, x_query: np.ndarray
) -> np.ndarray:
    if name == "linear":
        mapping = LinearRegression().fit(x_train.reshape(-1, 1), y_train)
        prediction = mapping.predict(x_query.reshape(-1, 1))
        return np.clip(prediction, y_train.min(), y_train.max())
    if name == "isotonic":
        mapping = IsotonicRegression(increasing=False, out_of_bounds="clip")
        return mapping.fit(x_train, y_train).predict(x_query)
    raise ValueError(name)


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    direct = (
        pl.read_parquet(inp / "measurements.parquet")
        .filter(pl.col("endpoint") == cfg["direct_endpoint"])
        .select("compound_id", pl.col("value").alias("y_true"), "conf_low", "conf_high")
    )
    primary = (
        pl.read_parquet(inp / "regression-labels.parquet")
        .filter(pl.col("head") == cfg["primary_head"])
        .select("compound_id", pl.col("value").alias("primary_value"))
    )
    scores = pl.read_parquet(inp / "scores.parquet")
    medium = pl.read_parquet(inp / "medium.parquet")
    feature_cols = sorted(c for c in scores.columns if c.startswith("score_")) + sorted(
        c for c in medium.columns if c.startswith("emb_")
    )
    features = scores.select(
        "compound_id", *[c for c in feature_cols if c.startswith("score_")]
    ).join(
        medium.select("compound_id", *[c for c in feature_cols if c.startswith("emb_")]),
        on="compound_id",
        validate="1:1",
    )
    if direct.height != 1493 or primary.height != 4376 or len(feature_cols) != 647:
        raise ValueError("Unexpected label or feature coverage")
    output = []
    map_rows = []
    for split in cfg["splits"]:
        assignments = pl.read_parquet(inp / f"{split}-assignments.parquet")
        membership = assignments.select("compound_id", "fold", "is_anchor")
        data = (
            membership.join(features, on="compound_id", validate="1:1")
            .join(direct, on="compound_id", how="left", validate="1:1")
            .join(primary, on="compound_id", how="left", validate="1:1")
            .sort("compound_id")
        )
        if data.filter(pl.col("y_true").is_not_null()).height != 1493:
            raise ValueError(f"Incomplete {split} direct labels")
        fold_ids = data["fold"].to_numpy()
        anchors = data["is_anchor"].to_numpy()
        y = data["y_true"].to_numpy()
        z = data["primary_value"].to_numpy()
        x = data.select(feature_cols).to_numpy().astype(np.float32)
        if not np.isfinite(x).all():
            raise ValueError(f"Nonfinite {split} features")
        for fold in range(5):
            allowed = (fold_ids != fold) | anchors
            truth_train = allowed & np.isfinite(y)
            proxy_train = allowed & np.isfinite(z) & ~np.isfinite(y)
            valid = (fold_ids == fold) & ~anchors & np.isfinite(y)
            paired = truth_train & np.isfinite(z)
            if sum(valid) < 100 or sum(proxy_train) < 1500 or sum(paired) < 500:
                raise ValueError(f"Insufficient {split}/fold{fold} coverage")
            variants = [("direct_only", None, 0.0)] + [
                (f"{mapping}_weight_{weight:g}", mapping, weight)
                for mapping in cfg["mappings"]
                for weight in cfg["pseudo_weights"]
            ]
            pseudo_values = {}
            for mapping in cfg["mappings"]:
                pseudo_values[mapping] = map_primary(mapping, z[paired], y[paired], z[proxy_train])
                observed_valid = valid & np.isfinite(z)
                mapped_valid = map_primary(mapping, z[paired], y[paired], z[observed_valid])
                map_rows.append(
                    {
                        "split_id": split,
                        "fold": fold,
                        "mapping": mapping,
                        "paired_train": int(sum(paired)),
                        "primary_only_train": int(sum(proxy_train)),
                        "valid_paired": int(sum(observed_valid)),
                        "oracle_primary_mae": float(
                            np.abs(y[observed_valid] - mapped_valid).mean()
                        ),
                        "oracle_primary_spearman": float(
                            spearmanr(y[observed_valid], mapped_valid).statistic
                        ),
                    }
                )
            for name, mapping, weight in variants:
                if mapping is None:
                    fit_x = x[truth_train]
                    fit_y = y[truth_train]
                    fit_weight = np.ones(sum(truth_train))
                else:
                    fit_x = np.concatenate([x[truth_train], x[proxy_train]])
                    fit_y = np.concatenate([y[truth_train], pseudo_values[mapping]])
                    fit_weight = np.concatenate(
                        [np.ones(sum(truth_train)), np.full(sum(proxy_train), weight)]
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
                model.fit(fit_x, fit_y, sample_weight=fit_weight)
                pred = model.predict(x[valid])
                if not np.isfinite(pred).all():
                    raise ValueError(f"Nonfinite {split}/{fold}/{name} prediction")
                output.append(
                    data.filter(pl.Series(valid))
                    .select("compound_id", "fold", "y_true", "conf_low", "conf_high")
                    .with_columns(
                        pl.lit(cfg["direct_endpoint"]).alias("endpoint"),
                        pl.lit(split).alias("split_id"),
                        pl.lit(name).alias("method"),
                        pl.Series("prediction", pred),
                    )
                )
                print(split, fold, name, flush=True)
    oof = pl.concat(output).sort("split_id", "method", "compound_id")
    expected = 2 * (1 + len(cfg["mappings"]) * len(cfg["pseudo_weights"])) * 1457
    if (
        oof.height != expected
        or oof.select("split_id", "method", "compound_id").is_duplicated().any()
    ):
        raise ValueError("Incomplete OOF output")
    metrics = []
    for (split, method), frame in oof.group_by("split_id", "method"):
        for fold in [-1, 0, 1, 2, 3, 4]:
            part = frame if fold == -1 else frame.filter(pl.col("fold") == fold)
            yt, lower, upper, pred = (
                part[c].to_numpy() for c in ["y_true", "conf_low", "conf_high", "prediction"]
            )
            metrics.append(
                {
                    "split_id": split,
                    "method": method,
                    "fold": fold,
                    "n": part.height,
                    "ST_RAE": st_rae(yt, lower, upper, pred),
                    "MAE": float(np.abs(yt - pred).mean()),
                    "Spearman": float(spearmanr(yt, pred).statistic),
                    "weak_MAE": float(np.abs(yt[yt < 4] - pred[yt < 4]).mean()),
                }
            )
    oof.write_parquet(out / "oof.parquet")
    pl.DataFrame(metrics).sort("split_id", "method", "fold").write_csv(out / "metrics.csv")
    pl.DataFrame(map_rows).sort("split_id", "mapping", "fold").write_csv(
        out / "mapping-diagnostics.csv"
    )
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "oof_rows": oof.height,
                "qualification": "Exploratory proxy-supervision screen; outer-fold exclusion",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
