"""Audit matched raw and PCA CheMeleon TabICL OOF predictions."""

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
value_columns = ["y_true", "conf_low", "conf_high", "is_neighbor"]


def read(name: str, prediction_name: str) -> pl.DataFrame:
    return (
        pl.read_parquet(inp / name)
        .sort(keys)
        .select(*keys, *value_columns, pl.col("prediction").alias(prediction_name))
    )


raw = read(cfg["raw"], "raw2048")
pca = read(cfg["pca"], "pca256")
assert raw.height == pca.height == 6384
assert raw.select(keys).n_unique() == pca.select(keys).n_unique() == 6384
for column in keys + value_columns:
    assert raw[column].equals(pca[column]), column

aligned = raw.join(pca.select(*keys, "pca256"), on=keys, validate="1:1").with_columns(
    pl.mean_horizontal("raw2048", "pca256").alias("raw2048_pca256_equal")
)
methods = ["raw2048", "pca256", "raw2048_pca256_equal"]
assert aligned.select(methods).select(pl.all().is_finite().all()).row(0) == (True,) * 3

spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)

rows = []
for (endpoint,), endpoint_table in aligned.group_by("endpoint", maintain_order=True):
    for fold in [-1, 0, 1, 2, 3, 4]:
        table = endpoint_table if fold == -1 else endpoint_table.filter(pl.col("fold") == fold)
        y = table["y_true"].to_numpy()
        for method in methods:
            prediction = table[method].to_numpy()
            rows.append(
                {
                    "endpoint": endpoint,
                    "fold": fold,
                    "n": table.height,
                    "method": method,
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

metrics = pl.DataFrame(rows).sort("endpoint", "fold", "method")
pooled = metrics.filter(pl.col("fold") == -1)
macro = pooled.group_by("method").agg(pl.col("ST_RAE", "MAE", "Spearman").mean()).sort("ST_RAE")
fold_deltas = (
    metrics.filter(pl.col("fold") >= 0)
    .pivot(on="method", index=["endpoint", "fold"], values="ST_RAE")
    .with_columns(
        (pl.col("raw2048") - pl.col("pca256")).alias("raw_minus_pca"),
        (pl.col("raw2048_pca256_equal") - pl.col("pca256")).alias("equal_minus_pca"),
    )
    .sort("endpoint", "fold")
)
correlations = aligned.group_by("endpoint").agg(
    pl.corr("raw2048", "pca256", method="pearson").alias("prediction_pearson"),
    pl.corr("raw2048", "pca256", method="spearman").alias("prediction_spearman"),
    pl.corr(
        pl.col("raw2048") - pl.col("y_true"),
        pl.col("pca256") - pl.col("y_true"),
        method="pearson",
    ).alias("residual_pearson"),
)

aligned.write_parquet(out / "aligned-oof.parquet")
metrics.write_csv(out / "metrics.csv")
macro.write_csv(out / "macro.csv")
fold_deltas.write_csv(out / "fold-deltas.csv")
correlations.write_csv(out / "correlations.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "rows": aligned.height,
            "methods": methods,
            "weight_search": False,
            "ensemble": "Equal mean of raw2048 and PCA256 predictions",
            "macro": macro.to_dicts(),
        },
        indent=2,
    )
)
print(macro)
print(pooled.sort("endpoint", "ST_RAE"))
print(correlations.sort("endpoint"))
