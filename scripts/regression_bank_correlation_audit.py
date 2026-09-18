"""Audit atomic regression-bank performance and pairwise prediction redundancy."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)

metric_rows = []
pair_rows = []
keys = ["compound_id", "endpoint", "fold"]
for target, pool in cfg["pools"].items():
    frames = []
    for name, file in pool["components"]:
        frame = (
            pl.read_parquet(inp / file)
            .filter(pl.col("endpoint") == pool["endpoint"])
            .sort(keys)
        )
        if frame.height == 0 or frame.select(keys).n_unique() != frame.height:
            raise ValueError(f"Invalid atomic OOF: {target}/{name}")
        frames.append((name, frame))
    reference = frames[0][1]
    y = reference["y_true"].to_numpy()
    low = reference["conf_low"].to_numpy()
    high = reference["conf_high"].to_numpy()
    predictions = {}
    for name, frame in frames:
        for column in keys + ["y_true", "conf_low", "conf_high"]:
            if not frame[column].equals(reference[column]):
                raise ValueError(f"Misaligned {target}/{name}/{column}")
        prediction = frame["prediction"].to_numpy()
        predictions[name] = prediction
        metric_rows.append(
            {
                "target": target,
                "method": name,
                "ST_RAE": float(
                    scorer.rae_soft_threshold_absolute_error(
                        y,
                        prediction,
                        y_true_lower=low,
                        y_true_upper=high,
                    )
                ),
                "MAE": float(np.abs(y - prediction).mean()),
                "Spearman_y": float(spearmanr(y, prediction).statistic),
            }
        )
    names = list(predictions)
    for left_index, left in enumerate(names):
        for right in names[left_index + 1 :]:
            left_prediction = predictions[left]
            right_prediction = predictions[right]
            pair_rows.append(
                {
                    "target": target,
                    "left": left,
                    "right": right,
                    "Pearson_prediction": float(
                        np.corrcoef(left_prediction, right_prediction)[0, 1]
                    ),
                    "Spearman_prediction": float(
                        spearmanr(left_prediction, right_prediction).statistic
                    ),
                    "Pearson_error": float(
                        np.corrcoef(y - left_prediction, y - right_prediction)[0, 1]
                    ),
                    "mean_absolute_difference": float(
                        np.abs(left_prediction - right_prediction).mean()
                    ),
                }
            )

metrics = pl.DataFrame(metric_rows).sort("target", "ST_RAE")
pairs = pl.DataFrame(pair_rows).sort(
    "target", "Spearman_prediction", descending=[False, True]
)
redundant = pairs.filter(pl.col("Spearman_prediction") >= 0.985)
metrics.write_csv(out / "atomic-metrics.csv")
pairs.write_csv(out / "pairwise-correlations.csv")
redundant.write_csv(out / "redundant-pairs-spearman-0985.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "redundancy_threshold": {"metric": "Spearman_prediction", "value": 0.985},
            "redundant_pairs": redundant.to_dicts(),
            "best_atomic_by_target": metrics.group_by(
                "target", maintain_order=True
            ).head(1).to_dicts(),
        },
        indent=2,
    )
    + "\n"
)
print(redundant, flush=True)
