"""Align nested auxiliary predictions to corrected outer validation compounds."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
base = pl.read_parquet(inp / cfg["base"]).filter(
    (pl.col("mode") == "parent_retained") & (pl.col("features") == "full")
)
true = pl.read_parquet(inp / cfg["primary"]).with_columns(
    pl.col("head").str.split("@").list.first().alias("enzyme")
)
aux = []
for item in cfg["aux"]:
    f = pl.concat([pl.read_parquet(inp / p) for p in item["predictions"]])
    for fold in item["excluded_folds"]:
        for head in [c for c in f.columns if c != "compound_id"]:
            aux.append(
                f.select("compound_id", pl.col(head).alias("predicted_primary")).with_columns(
                    pl.lit(fold).alias("fold"), pl.lit(head.split("@")[0]).alias("enzyme")
                )
            )
pred = (
    pl.concat(aux)
    .group_by("compound_id", "fold", "enzyme")
    .agg(pl.col("predicted_primary").mean(), pl.len().alias("nmodels"))
)
assert (pred["nmodels"] == 4).all()
f = (
    base.with_columns(pl.col("endpoint").str.split("_").list.first().alias("enzyme"))
    .join(pred, on=["compound_id", "fold", "enzyme"], validate="m:1")
    .join(
        true.select("compound_id", "enzyme", pl.col("value").alias("observed_primary")),
        on=["compound_id", "enzyme"],
        how="left",
        validate="m:1",
    )
    .filter(pl.col("observed_primary").is_not_null())
    .with_columns(
        (pl.col("predicted_primary") - pl.col("observed_primary")).alias("primary_residual"),
        (pl.col("prediction") - pl.col("y_true")).alias("pic50_residual"),
    )
)
rows = []
for ep, g in f.group_by("endpoint"):
    for subset in ["all", "neighbors"]:
        h = g if subset == "all" else g.filter(pl.col("is_neighbor"))
        rows.append(
            dict(
                endpoint=ep[0],
                subset=subset,
                n=h.height,
                primary_RMSE=float(np.sqrt((h["primary_residual"].to_numpy() ** 2).mean())),
                primary_Spearman=float(
                    spearmanr(h["predicted_primary"], h["observed_primary"]).statistic
                ),
                residual_Spearman=float(
                    spearmanr(h["primary_residual"], h["pic50_residual"]).statistic
                ),
            )
        )
f.write_parquet(out / "aligned-residuals.parquet")
pl.DataFrame(rows).write_csv(out / "metrics.csv")
write_json(
    out / "summary.json",
    dict(
        config=cfg,
        scope="Descriptive residual association; no hidden primary measurements.",
    ),
)
print(pl.DataFrame(rows))
