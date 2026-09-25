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
    split = (
        pl.read_parquet(inp / cfg["assignments"])
        .filter(pl.col("strategy") == "umap")
        .select("compound_id", "fold")
    )
    regression = pl.read_parquet(inp / cfg["regression_labels"]).filter(
        pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
    )
    classification = pl.read_parquet(inp / cfg["classification_labels"]).filter(
        pl.col("strategy") == "umap"
    )
    labels = [
        ("regression", regression.join(split, on="compound_id", validate="m:1")),
        ("classification", classification.select("compound_id", "endpoint", "value", "fold")),
    ]
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
                endpoint_predictions = np.empty(len(y), dtype=np.float64)
                for fold in range(5):
                    train = folds != fold
                    if track == "regression":
                        model = lgb.LGBMRegressor(
                            objective="regression_l1", **cfg["lightgbm"], random_state=42 + fold
                        )
                        model.fit(x[train], y[train])
                        values = model.predict(x[~train])
                    else:
                        model = lgb.LGBMClassifier(
                            objective="binary", **cfg["lightgbm"], random_state=42 + fold
                        )
                        model.fit(x[train], y[train].astype(int))
                        values = model.predict_proba(x[~train])[:, 1]
                    endpoint_predictions[~train] = values
                    print(
                        f"depth={depth} track={track} endpoint={endpoint} fold={fold}",
                        flush=True,
                    )
                assert np.isfinite(endpoint_predictions).all()
                predictions.append(
                    table.select("compound_id", "endpoint", "fold")
                    .with_columns(
                        pl.lit(track).alias("track"),
                        pl.lit(depth).alias("depth"),
                        pl.lit(cfg["split_id"]).alias("split_id"),
                        pl.Series("y_true", y),
                        pl.Series("prediction", endpoint_predictions),
                    )
                )
                row = {"track": track, "endpoint": endpoint, "depth": depth, "n": len(y)}
                if track == "regression":
                    row.update(
                        ST_RAE=float(
                            scorer.rae_soft_threshold_absolute_error(
                                y,
                                endpoint_predictions,
                                y_true_lower=table["conf_low"].to_numpy(),
                                y_true_upper=table["conf_high"].to_numpy(),
                            )
                        ),
                        MAE=float(np.abs(y - endpoint_predictions).mean()),
                        Spearman=float(spearmanr(y, endpoint_predictions).statistic),
                    )
                else:
                    row.update(
                        AUROC=float(roc_auc_score(y, endpoint_predictions)),
                        AP=float(average_precision_score(y, endpoint_predictions)),
                    )
                metrics.append(row)
    pl.concat(predictions, how="diagonal").write_parquet(out / "oof.parquet")
    pl.DataFrame(metrics).write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(
        json.dumps({"config": cfg, "metrics": metrics}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
