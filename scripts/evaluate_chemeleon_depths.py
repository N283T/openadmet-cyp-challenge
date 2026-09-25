"""Matched UMAP CV comparison of frozen CheMeleon depths on both CYP tracks."""

import importlib.util
import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score


def scorer_module(path):
    spec = importlib.util.spec_from_file_location("official_scoring", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    scorer = scorer_module(inp / cfg["scorer"])
    assignments = pl.read_parquet(inp / cfg["assignments"]).filter(
        pl.col("strategy") == cfg.get("strategy", "umap")
    )
    if "is_anchor" not in assignments.columns:
        assignments = assignments.with_columns(pl.lit(False).alias("is_anchor"))
    split = assignments.select("compound_id", "fold", "is_anchor")
    labels = []
    if "regression" in cfg.get("tracks", ["regression", "classification"]):
        regression = pl.read_parquet(inp / cfg["regression_labels"]).filter(
            pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
        )
        labels.append(("regression", regression.join(split, on="compound_id", validate="m:1")))
    if "classification" in cfg.get("tracks", ["regression", "classification"]):
        classification = pl.read_parquet(inp / cfg["classification_labels"]).filter(
            pl.col("strategy") == "umap"
        )
        labels.append(
            (
                "classification",
                classification.select("compound_id", "endpoint", "value", "fold").with_columns(
                    pl.lit(False).alias("is_anchor")
                ),
            )
        )
    predictions = []
    metrics = []
    for depth in cfg["depths"]:
        feature = pl.read_parquet(inp / cfg["features"][str(depth)])
        columns = [name for name in feature.columns if name.startswith("emb_")]
        assert len(columns) == cfg["dimensions"]
        assert feature["compound_id"].n_unique() == feature.height
        for track, frame in labels:
            for endpoint in sorted(frame["endpoint"].unique()):
                table = (
                    frame.filter(pl.col("endpoint") == endpoint)
                    .join(feature, on="compound_id", validate="m:1")
                    .sort("compound_id")
                )
                assert table["fold"].null_count() == 0
                x = table.select(columns).to_numpy().astype(np.float32)
                y = table["value"].to_numpy()
                folds = table["fold"].to_numpy()
                assert np.isfinite(x).all() and np.isfinite(y).all()
                anchors = table["is_anchor"].to_numpy()
                endpoint_predictions = np.full(len(y), np.nan, dtype=np.float64)
                for fold in range(5):
                    validation = (folds == fold) & ~anchors
                    train = ~validation
                    if track == "regression":
                        model = lgb.LGBMRegressor(
                            objective="regression_l1", **cfg["lightgbm"], random_state=42 + fold
                        )
                        model.fit(x[train], y[train])
                        values = model.predict(x[validation])
                    else:
                        model = lgb.LGBMClassifier(
                            objective="binary", **cfg["lightgbm"], random_state=42 + fold
                        )
                        model.fit(x[train], y[train].astype(int))
                        values = model.predict_proba(x[validation])[:, 1]
                    endpoint_predictions[validation] = values
                    print(
                        f"depth={depth} track={track} endpoint={endpoint} fold={fold}",
                        flush=True,
                    )
                evaluated = ~anchors
                assert np.isfinite(endpoint_predictions[evaluated]).all()
                assert np.isnan(endpoint_predictions[anchors]).all()
                observed = table.filter(pl.Series(evaluated))
                observed_y = y[evaluated]
                observed_prediction = endpoint_predictions[evaluated]
                predictions.append(
                    observed.select("compound_id", "endpoint", "fold").with_columns(
                        pl.lit(track).alias("track"),
                        pl.lit(depth).alias("depth"),
                        pl.lit(cfg["split_id"]).alias("split_id"),
                        pl.Series("y_true", observed_y),
                        pl.Series("prediction", observed_prediction),
                    )
                )
                row = {
                    "track": track,
                    "endpoint": endpoint,
                    "depth": depth,
                    "n": len(observed_y),
                }
                if track == "regression":
                    row.update(
                        ST_RAE=float(
                            scorer.rae_soft_threshold_absolute_error(
                                observed_y,
                                observed_prediction,
                                y_true_lower=observed["conf_low"].to_numpy(),
                                y_true_upper=observed["conf_high"].to_numpy(),
                            )
                        ),
                        MAE=float(np.abs(observed_y - observed_prediction).mean()),
                        Spearman=float(spearmanr(observed_y, observed_prediction).statistic),
                    )
                else:
                    row.update(
                        AUROC=float(roc_auc_score(observed_y, observed_prediction)),
                        AP=float(average_precision_score(observed_y, observed_prediction)),
                    )
                metrics.append(row)
    pl.concat(predictions, how="diagonal").write_parquet(out / "oof.parquet")
    pl.DataFrame(metrics).write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(
        json.dumps({"config": cfg, "metrics": metrics}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
