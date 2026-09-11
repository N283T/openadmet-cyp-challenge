"""Matched OOF prediction/residual correlations and prespecified equal blends."""

import importlib.util
import itertools
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import psycopg
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
frame = pl.concat(pl.read_parquet(inp / p) for p in cfg["oofs"])
keys = ["compound_id", "endpoint", "fold"]
methods = sorted(frame["method"].unique())
assert len(methods) == 3
base = frame.filter(pl.col("method") == methods[0]).sort(keys)
assert base.height == 6384 and not base.select(keys).is_duplicated().any()
aligned = base.select(keys + ["split_id", "y_true", "conf_low", "conf_high", "is_neighbor"])
for method in methods:
    h = frame.filter(pl.col("method") == method).sort(keys)
    for col in keys + ["split_id", "y_true", "conf_low", "conf_high", "is_neighbor"]:
        assert h[col].equals(aligned[col]), (method, col)
    assert h["prediction"].is_finite().all()
    aligned = aligned.with_columns(h["prediction"].alias(method))
# Controls have the same labels and folds but differently named split identifiers.
for method, path in cfg["controls"].items():
    h = pl.read_parquet(inp / path).sort(keys)
    for col in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
        assert h[col].equals(base[col]), (method, col)
    aligned = aligned.with_columns(h["prediction"].alias(method))
blends = {}
for size in [2, 3]:
    for combo in itertools.combinations(methods, size):
        name = "mean__" + "__".join(combo)
        aligned = aligned.with_columns(pl.mean_horizontal([pl.col(c) for c in combo]).alias(name))
        blends[name] = combo
aligned.write_parquet(out / "aligned-oof.parquet")
correlations, metrics = [], []
for (ep,), group in aligned.group_by("endpoint"):
    for subset in ["all", "neighbors"]:
        sub = group if subset == "all" else group.filter(pl.col("is_neighbor"))
        for fold in [-1, 0, 1, 2, 3, 4]:
            h = sub if fold == -1 else sub.filter(pl.col("fold") == fold)
            y = h["y_true"].to_numpy()
            for a, b in itertools.combinations(methods, 2):
                pa, pb = h[a].to_numpy(), h[b].to_numpy()
                correlations.append(
                    dict(
                        endpoint=ep,
                        subset=subset,
                        fold=fold,
                        n=h.height,
                        a=a,
                        b=b,
                        prediction_pearson=float(np.corrcoef(pa, pb)[0, 1]),
                        prediction_spearman=float(spearmanr(pa, pb).statistic),
                        residual_pearson=float(np.corrcoef(pa - y, pb - y)[0, 1]),
                    )
                )
            for method in methods + list(blends) + list(cfg["controls"]):
                p = h[method].to_numpy()
                metrics.append(
                    dict(
                        endpoint=ep,
                        subset=subset,
                        fold=fold,
                        n=h.height,
                        method=method,
                        ST_RAE=float(
                            scorer.rae_soft_threshold_absolute_error(
                                y,
                                p,
                                y_true_lower=h["conf_low"].to_numpy(),
                                y_true_upper=h["conf_high"].to_numpy(),
                            )
                        ),
                        MAE=float(np.abs(y - p).mean()),
                        Spearman=float(spearmanr(y, p).statistic),
                    )
                )
pl.DataFrame(correlations).sort("subset", "endpoint", "fold", "a", "b").write_csv(
    out / "correlations.csv"
)
m = pl.DataFrame(metrics)
m.write_csv(out / "metrics.csv")
macro = (
    m.filter(pl.col("fold") == -1)
    .group_by("method", "subset")
    .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
    .sort("subset", "ST_RAE")
)
macro.write_csv(out / "macro.csv")
with psycopg.connect(**cfg["database"]) as conn:
    for r in metrics:
        for name in ["ST_RAE", "MAE", "Spearman"]:
            conn.execute(
                "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
                (
                    Path(os.environ["CYP_RUN_DIR"]).name,
                    r["endpoint"],
                    "ensemble_components:" + r["subset"],
                    str(r["fold"]),
                    r["method"] + ":" + name,
                    r[name],
                ),
            )
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            matched_oof_rows=base.height,
            blends=blends,
            weight_selection="none; equal pairs and equal triple",
        ),
        indent=2,
    )
)
print(macro)
