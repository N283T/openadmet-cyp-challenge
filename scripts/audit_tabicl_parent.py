"""Audit TabICL parent-retained baselines against matched LightGBM predictions."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import psycopg
from scipy.stats import spearmanr

from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
oof = pl.concat([pl.read_parquet(inp / p) for p in cfg["oof"]]).with_columns(
    (pl.lit("tabicl_") + pl.col("features")).alias("method")
)
base = (
    pl.read_parquet(inp / cfg["baseline"])
    .filter((pl.col("features") == "full") & (pl.col("mode") == "parent_retained"))
    .with_columns(pl.lit("lightgbm_full").alias("method"))
)
assert oof.height == 2 * base.height
oof = pl.concat([oof, base], how="vertical_relaxed")
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
official = importlib.util.module_from_spec(spec)
spec.loader.exec_module(official)
rows = []
for (ep, method), f in oof.group_by("endpoint", "method"):
    ref = base.filter(pl.col("endpoint") == ep).select(
        "compound_id", "fold", pl.col("y_true").alias("reference_y")
    )
    assert f.height == ref.height
    f = f.join(ref, on=["compound_id", "fold"], validate="1:1", maintain_order="left")
    assert np.allclose(f["y_true"], f["reference_y"])
    for subset in ["all_nonanchors", "neighbors"]:
        g = f if subset == "all_nonanchors" else f.filter(pl.col("is_neighbor"))
        for fold in [-1, 0, 1, 2, 3, 4]:
            h = g if fold == -1 else g.filter(pl.col("fold") == fold)
            y, p, lo, hi = [
                h[c].to_numpy() for c in ["y_true", "prediction", "conf_low", "conf_high"]
            ]
            assert np.isfinite(p).all()
            rows.append(
                dict(
                    endpoint=ep,
                    method=method,
                    subset=subset,
                    fold=fold,
                    n=h.height,
                    ST_RAE=float(
                        official.rae_soft_threshold_absolute_error(
                            y, p, y_true_lower=lo, y_true_upper=hi
                        )
                    ),
                    MAE=float(np.abs(y - p).mean()),
                    Spearman=float(spearmanr(y, p).statistic),
                )
            )
m = pl.DataFrame(rows)
with psycopg.connect(**cfg["database"]) as conn:
    for r in rows:
        for metric in ["ST_RAE", "MAE", "Spearman"]:
            conn.execute(
                "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
                (
                    Path(os.environ["CYP_RUN_DIR"]).name,
                    r["endpoint"],
                    "parent_tabicl:" + r["subset"],
                    str(r["fold"]),
                    r["method"] + ":" + metric,
                    r[metric],
                ),
            )
m.write_csv(out / "metrics.csv")
macro = (
    m.filter(pl.col("fold") == -1)
    .group_by("method", "subset")
    .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
    .sort("subset", "ST_RAE")
)
macro.write_csv(out / "macro.csv")
oof.write_parquet(out / "oof.parquet")
write_json(
    out / "summary.json",
    dict(
        config=cfg,
        results=macro.to_dicts(),
        matched_key_label_checks=True,
        scope="Parent CV, fixed means, matched validation; compact256+4 and full+Mordred.",
    ),
)
print(macro)
