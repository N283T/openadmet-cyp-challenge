"""Keyed matched-learner audit using the official ST-RAE scorer."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import psycopg
from scipy.stats import spearmanr

from cyp_challenge.calibration_alignment import aligned_matrix
from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
labels = pl.read_parquet(inp / cfg["labels"])
oof = pl.concat([pl.read_parquet(inp / p) for p in cfg["oof"]])
test = pl.concat([pl.read_parquet(inp / p) for p in cfg["test"]])
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
assert oof.select("compound_id", "endpoint", "method").n_unique() == oof.height == 6525 * 3
reference = pl.read_parquet(inp / cfg["distribution"])
rows = []
for ep in sorted(labels["endpoint"].unique()):
    f = labels.filter(pl.col("endpoint") == ep).sort("compound_id")
    y, lo, hi = [f[c].to_numpy() for c in ["y_true", "conf_low", "conf_high"]]
    predictions = {"tabpfn_mean": f["prediction"].to_numpy()}
    for method in ["baseline:median", "baseline:q0.55", "blend3:q0.55"]:
        part = reference.filter(
            (pl.col("endpoint") == ep) & (pl.col("method") == method) & (pl.col("repeat") == -1)
        )
        predictions[method] = aligned_matrix(f, part, ["prediction"])[:, 0]
    for method in sorted(oof["method"].unique()):
        p = oof.filter((pl.col("endpoint") == ep) & (pl.col("method") == method))
        predictions[method] = aligned_matrix(f, p, ["prediction"])[:, 0]
    for method, p in predictions.items():
        for fold in [-1, 0, 1, 2, 3, 4]:
            ix = np.ones(len(y), dtype=bool) if fold == -1 else f["fold"].to_numpy() == fold
            rows.append(
                dict(
                    endpoint=ep,
                    method=method,
                    fold=fold,
                    ST_RAE=float(
                        scorer.rae_soft_threshold_absolute_error(
                            y[ix], p[ix], y_true_lower=lo[ix], y_true_upper=hi[ix]
                        )
                    ),
                    MAE=float(np.abs(y[ix] - p[ix]).mean()),
                    Spearman=float(spearmanr(y[ix], p[ix]).statistic),
                )
            )
metrics = pl.DataFrame(rows)
summary = (
    metrics.filter(pl.col("fold") == -1)
    .group_by("method")
    .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
    .sort("ST_RAE")
)
with psycopg.connect(**cfg["database"]) as conn:
    for r in metrics.iter_rows(named=True):
        for metric in ["ST_RAE", "MAE", "Spearman"]:
            conn.execute(
                "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
                (
                    Path(os.environ["CYP_RUN_DIR"]).name,
                    r["endpoint"],
                    "umap:matched-tabicl",
                    str(r["fold"]),
                    r["method"] + ":" + metric,
                    r[metric],
                ),
            )
metrics.write_csv(out / "metrics.csv")
summary.write_csv(out / "summary.csv")
oof.write_parquet(out / "oof.parquet")
test = test.group_by("compound_id", "endpoint", "method").agg(
    pl.col("prediction").mean(), pl.col("fold").n_unique().alias("nfold")
)
assert (test["nfold"] == 5).all()
test.write_parquet(out / "test.parquet")
test.group_by("endpoint", "method").agg(
    pl.col("prediction").mean().alias("mean"),
    pl.col("prediction").std().alias("sd"),
    (pl.col("prediction") < 4).sum().alias("below4"),
).write_csv(out / "test-summary.csv")
write_json(
    out / "summary.json",
    dict(
        results=summary.to_dicts(),
        config=cfg,
        scope="Fixed UMAP comparison; official wrapper does not pin TabICL version.",
    ),
)
print(summary)
