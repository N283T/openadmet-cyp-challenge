"""Pooled auxiliary OOF evaluation; no model or hyperparameter selection."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

from cyp_challenge.storage import write_json


def main():
    inputs, output = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    predictions = pl.concat(
        [pl.read_parquet(p) for p in sorted(inputs.glob("*-outer_predictions.parquet"))]
    )
    targets = pl.read_parquet(next(inputs.glob("*-targets.parquet"))).rename({"head": "endpoint"})
    summaries = [json.loads(p.read_text()) for p in sorted(inputs.glob("*-summary.json"))]
    assert len(summaries) == len(config["expected_outer_folds"])
    assert {s["config"]["outer_fold"] for s in summaries} == set(config["expected_outer_folds"])
    assert predictions.select("compound_id", "endpoint").n_unique() == predictions.height
    assert predictions.height == targets.height
    joined = targets.join(predictions, on=["compound_id", "endpoint"], how="left", validate="1:1")
    assert joined["prediction"].null_count() == 0
    assert np.allclose(
        joined["value"].to_numpy(),
        joined["y_true"].to_numpy(),
        rtol=config["label_rtol"],
        atol=config["label_atol"],
    )
    rows = []
    for (endpoint,), frame in joined.group_by("endpoint"):
        y, p = frame["y_true"].to_numpy(), frame["prediction"].to_numpy()
        rows.append(
            {
                "endpoint": endpoint,
                "n": len(y),
                "rmse": float(np.sqrt(((y - p) ** 2).mean())),
                "mae": float(np.abs(y - p).mean()),
                "spearman": float(spearmanr(y, p).statistic),
            }
        )
    scores = pl.DataFrame(rows).sort("endpoint")
    scores.write_csv(output / "pooled_metrics.csv")
    predictions.write_parquet(output / "oof.parquet")
    write_json(
        output / "audit.json",
        {
            "passed": True,
            "prediction_rows": predictions.height,
            "metrics": scores.to_dicts(),
            "folds": [
                {
                    "fold": s["config"]["outer_fold"],
                    "best_epoch": s["best_epoch"],
                    "epochs": s["epochs"],
                    "counts": s["counts"],
                }
                for s in summaries
            ],
        },
    )
    print(scores, flush=True)


if __name__ == "__main__":
    main()
