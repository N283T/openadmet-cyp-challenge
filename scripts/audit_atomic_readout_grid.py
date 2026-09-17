"""Score a target-specific grid of atomic OOF readouts without blending them."""

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

rows = []
aligned = None
keys = ["compound_id", "endpoint", "fold"]
label_columns = ["y_true", "conf_low", "conf_high", "is_neighbor"]
for variant in cfg["variants"]:
    frame = pl.read_parquet(inp / variant["file"]).sort(keys)
    expected_endpoints = cfg["endpoints"] if "endpoints" in cfg else [cfg["endpoint"]]
    if "endpoint" in cfg:
        frame = frame.filter(pl.col("endpoint") == cfg["endpoint"])
    assert sorted(frame["endpoint"].unique()) == sorted(expected_endpoints)
    assert frame.select(keys).n_unique() == frame.height
    assert np.isfinite(frame["prediction"].to_numpy()).all()
    if aligned is None:
        aligned = frame.select(*keys, *label_columns)
    else:
        for column in keys + label_columns:
            assert frame[column].equals(aligned[column]), (variant["name"], column)
    aligned = aligned.with_columns(frame["prediction"].alias(variant["name"]))
    for endpoint in expected_endpoints:
        endpoint_frame = frame.filter(pl.col("endpoint") == endpoint)
        for subset in ["all", "neighbors"]:
            selected = (
                endpoint_frame
                if subset == "all"
                else endpoint_frame.filter(pl.col("is_neighbor"))
            )
            for fold in [-1, 0, 1, 2, 3, 4]:
                table = selected if fold == -1 else selected.filter(pl.col("fold") == fold)
                y = table["y_true"].to_numpy()
                prediction = table["prediction"].to_numpy()
                rows.append(
                    {
                        "endpoint": endpoint,
                        "method": variant["name"],
                        "subset": subset,
                        "fold": fold,
                        "n": table.height,
                        "ST_RAE": float(
                            scorer.rae_soft_threshold_absolute_error(
                                y,
                                prediction,
                                y_true_lower=table["conf_low"].to_numpy(),
                                y_true_upper=table["conf_high"].to_numpy(),
                            )
                        ),
                        "MAE": float(np.abs(y - prediction).mean()),
                        "Spearman": float(spearmanr(y, prediction).statistic),
                    }
                )

metrics = pl.DataFrame(rows).sort("subset", "fold", "ST_RAE")
pooled = metrics.filter(pl.col("fold") == -1).sort("subset", "ST_RAE")
macro = (
    pooled.group_by("method", "subset")
    .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
    .sort("subset", "ST_RAE")
)
correlations = []
names = [variant["name"] for variant in cfg["variants"]]
for left_index, left in enumerate(names):
    for right in names[left_index + 1 :]:
        for endpoint in expected_endpoints:
            endpoint_frame = aligned.filter(pl.col("endpoint") == endpoint)
            correlations.append(
                {
                    "endpoint": endpoint,
                    "left": left,
                    "right": right,
                    "prediction_pearson": float(
                        endpoint_frame.select(pl.corr(left, right)).item()
                    ),
                    "residual_pearson": float(
                        endpoint_frame.select(
                            pl.corr(
                                pl.col(left) - pl.col("y_true"),
                                pl.col(right) - pl.col("y_true"),
                            )
                        ).item()
                    ),
                }
            )

metrics.write_csv(out / "metrics.csv")
pooled.write_csv(out / "pooled.csv")
macro.write_csv(out / "macro.csv")
pl.DataFrame(correlations).write_csv(out / "correlations.csv")
aligned.write_parquet(out / "aligned-oof.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "rows": aligned.height,
            "weight_search": False,
            "pooled": pooled.to_dicts(),
            "macro": macro.to_dicts(),
        },
        indent=2,
    )
    + "\n"
)
print(pooled, flush=True)
print(macro, flush=True)
