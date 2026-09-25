"""Prespecified equal OOF blends of frozen CheMeleon depth readouts."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
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
    source = pl.read_parquet(inp / cfg["oof"])
    assert not source.select("track", "endpoint", "depth", "compound_id").is_duplicated().any()
    intervals = pl.read_parquet(inp / cfg["measurements"]).select(
        "compound_id", "endpoint", "conf_low", "conf_high"
    )
    keys = ["track", "endpoint", "compound_id", "fold", "split_id", "y_true"]
    wide = source.pivot(on="depth", index=keys, values="prediction")
    assert wide.height * 3 == source.height
    assert wide.select("1", "3", "6").null_count().sum_horizontal().sum() == 0
    methods = {
        "depth_1": pl.col("1"),
        "depth_3": pl.col("3"),
        "depth_6": pl.col("6"),
        "depth_1_6_equal": (pl.col("1") + pl.col("6")) / 2,
        "depth_1_3_6_equal": (pl.col("1") + pl.col("3") + pl.col("6")) / 3,
    }
    predictions = pl.concat(
        [
            wide.select(
                *keys,
                pl.lit(name).alias("method"),
                expression.alias("prediction"),
            )
            for name, expression in methods.items()
        ]
    )
    assert predictions["prediction"].is_finite().all()
    predictions.write_parquet(out / "oof-blends.parquet")
    metrics = []
    for (track, endpoint, method), frame in predictions.group_by("track", "endpoint", "method"):
        frame = frame.sort("compound_id")
        y = frame["y_true"].to_numpy()
        pred = frame["prediction"].to_numpy()
        row = {"track": track, "endpoint": endpoint, "method": method, "n": len(y)}
        if track == "regression":
            joined = frame.join(
                intervals, on=["compound_id", "endpoint"], validate="1:1", maintain_order="left"
            )
            assert joined.height == frame.height
            row.update(
                ST_RAE=float(
                    scorer.rae_soft_threshold_absolute_error(
                        y,
                        pred,
                        y_true_lower=joined["conf_low"].to_numpy(),
                        y_true_upper=joined["conf_high"].to_numpy(),
                    )
                ),
                MAE=float(np.abs(y - pred).mean()),
            )
        else:
            row.update(
                AUROC=float(roc_auc_score(y, pred)),
                AP=float(average_precision_score(y, pred)),
            )
        metrics.append(row)
    pl.DataFrame(metrics).sort("track", "endpoint", "method").write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(
        json.dumps({"config": cfg, "metrics": metrics}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
