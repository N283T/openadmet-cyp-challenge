"""Matched fixed three-versus-four-model comparison after CatBoost readout."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import psycopg
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
keys = ["compound_id", "endpoint", "fold"]
base = None
for name, filename in cfg["models"].items():
    frame = pl.read_parquet(inp / filename).sort(keys)
    assert frame.height == 6384 and frame.select(keys).unique().height == 6384
    if base is None:
        base = frame.select(keys + ["y_true", "conf_low", "conf_high", "is_neighbor", "split_id"])
    for column in keys + ["y_true", "conf_low", "conf_high", "is_neighbor", "split_id"]:
        assert base[column].equals(frame[column]), column
    base = base.with_columns(frame["prediction"].alias(name))
base = base.with_columns(
    pl.mean_horizontal("caruana", "frozen_pfn", "ft_lgbm").alias("old_three_equal"),
    pl.mean_horizontal("caruana", "frozen_pfn", "ft_lgbm", "ft_catboost").alias("four_equal"),
)
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
methods = list(cfg["models"]) + ["old_three_equal", "four_equal"]
rows = []
for (endpoint,), group in base.group_by("endpoint"):
    for subset in ["all", "neighbors"]:
        selected = group if subset == "all" else group.filter(pl.col("is_neighbor"))
        for fold in [-1, 0, 1, 2, 3, 4]:
            data = selected if fold == -1 else selected.filter(pl.col("fold") == fold)
            y = data["y_true"].to_numpy()
            for method in methods:
                p = data[method].to_numpy()
                rows.append(
                    dict(
                        endpoint=endpoint,
                        subset=subset,
                        fold=fold,
                        method=method,
                        ST_RAE=float(
                            scorer.rae_soft_threshold_absolute_error(
                                y,
                                p,
                                y_true_lower=data["conf_low"].to_numpy(),
                                y_true_upper=data["conf_high"].to_numpy(),
                            )
                        ),
                        MAE=float(np.abs(y - p).mean()),
                        Spearman=float(spearmanr(y, p).statistic),
                    )
                )
metrics = pl.DataFrame(rows)
metrics.write_csv(out / "metrics.csv")
macro = (
    metrics.filter(pl.col("fold") == -1)
    .group_by("method", "subset")
    .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
    .sort("subset", "ST_RAE")
)
macro.write_csv(out / "macro.csv")
base.write_parquet(out / "aligned-oof.parquet")
base.select(keys + ["y_true", "conf_low", "conf_high", "is_neighbor", "split_id"]).with_columns(
    base["four_equal"].alias("prediction"),
    pl.lit("caruana_jeremy_pfn_ftlgbm_ftcatboost_equal").alias("method"),
).write_parquet(out / "external-oof.parquet")
with psycopg.connect(**cfg["database"]) as conn:
    for row in rows:
        for metric in ["ST_RAE", "MAE", "Spearman"]:
            conn.execute(
                "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
                (
                    Path(os.environ["CYP_RUN_DIR"]).name,
                    row["endpoint"],
                    "catboost_fixed:" + row["subset"],
                    str(row["fold"]),
                    row["method"] + ":" + metric,
                    row[metric],
                ),
            )
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            matched_rows=6384,
            weight_search=False,
            weights_three=[1 / 3] * 3,
            weights_four=[0.25] * 4,
        ),
        indent=2,
    )
)
print(macro, flush=True)
