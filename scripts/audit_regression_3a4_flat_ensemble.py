"""Compare flat CYP3A4 equal-member ensembles with the legacy nested recipe."""

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
keys = ["compound_id", "endpoint", "fold"]
truth = ["y_true", "conf_low", "conf_high", "is_neighbor"]


def component(name: str) -> pl.DataFrame:
    return (
        pl.read_parquet(inp / cfg[name])
        .filter(pl.col("endpoint") == cfg["endpoint"])
        .sort(keys)
        .select(*keys, pl.col("prediction").alias(name))
    )


top2_source = pl.read_parquet(inp / cfg["caruana_top2"])
current = (
    top2_source.filter(
        (pl.col("strategy") == "crossfold_topk")
        & (pl.col("weighting") == "weighted")
        & (pl.col("k") == 2)
    )
    .sort(keys)
    .select(*keys, *truth, pl.col("five_component_prediction").alias("current_nested5"))
)
equal_top2 = (
    top2_source.filter(
        (pl.col("strategy") == "crossfold_topk")
        & (pl.col("weighting") == "equal")
        & (pl.col("k") == 2)
    )
    .sort(keys)
    .select(*keys, pl.col("caruana_prediction").alias("caruana_top2_equal"))
)
table = current.join(equal_top2, on=keys, validate="1:1")
for name in ["frozen", "scores", "finetuned", "pca", "raw"]:
    table = table.join(component(name), on=keys, validate="1:1")
assert table.height == 2300 and table.select(keys).n_unique() == 2300

other_members = ["frozen", "scores", "finetuned", "pca"]
other_sum = pl.sum_horizontal(*other_members)
table = table.with_columns(
    (other_sum / 4).alias("flat_no_caruana4"),
    ((other_sum + pl.col("raw")) / 5).alias("flat_no_caruana_raw5"),
    ((2 * pl.col("caruana_top2_equal") + other_sum) / 6).alias("flat_caruana2_plus4"),
    ((2 * pl.col("caruana_top2_equal") + other_sum + pl.col("raw")) / 7).alias(
        "flat_caruana2_plus4_raw"
    ),
)
methods = [
    "current_nested5",
    "flat_no_caruana4",
    "flat_no_caruana_raw5",
    "flat_caruana2_plus4",
    "flat_caruana2_plus4_raw",
]
assert table.select(methods).select(pl.all().is_finite().all()).row(0) == (True,) * len(methods)

spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
rows = []
for fold in [-1, 0, 1, 2, 3, 4]:
    frame = table if fold == -1 else table.filter(pl.col("fold") == fold)
    y = frame["y_true"].to_numpy()
    for method in methods:
        prediction = frame[method].to_numpy()
        rows.append(
            {
                "fold": fold,
                "n": frame.height,
                "method": method,
                "ST_RAE": float(
                    scorer.rae_soft_threshold_absolute_error(
                        y,
                        prediction,
                        y_true_lower=frame["conf_low"].to_numpy(),
                        y_true_upper=frame["conf_high"].to_numpy(),
                    )
                ),
                "MAE": float(np.abs(y - prediction).mean()),
                "Spearman": float(spearmanr(y, prediction).statistic),
            }
        )

metrics = pl.DataFrame(rows).sort("fold", "ST_RAE")
pooled = metrics.filter(pl.col("fold") == -1).sort("ST_RAE")
fold_scores = metrics.filter(pl.col("fold") >= 0).pivot(
    on="method", index="fold", values="ST_RAE"
)
for method in methods[1:]:
    fold_scores = fold_scores.with_columns(
        (pl.col(method) - pl.col("current_nested5")).alias(f"{method}_minus_current")
    )

table.write_parquet(out / "aligned-oof.parquet")
metrics.write_csv(out / "metrics.csv")
pooled.write_csv(out / "pooled.csv")
fold_scores.write_csv(out / "fold-scores.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "rows": table.height,
            "methods": methods,
            "weight_search": False,
            "pooled": pooled.to_dicts(),
        },
        indent=2,
    )
)
print(pooled)
print(fold_scores)
