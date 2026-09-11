"""Verify complete updated foundation-model OOF and test predictions."""

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
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
official = importlib.util.module_from_spec(spec)
spec.loader.exec_module(official)
keys = ["compound_id", "endpoint", "fold"]
base = (
    pl.read_parquet(inp / cfg["lgbm"])
    .filter((pl.col("mode") == "parent_retained") & (pl.col("features") == "full"))
    .sort(keys)
)
frames = {"lgbm_updated": base}
frames["tabicl_previous"] = (
    pl.read_parquet(inp / cfg["previous_tabicl"])
    .filter(pl.col("method") == "tabicl_full")
    .sort(keys)
)
mapping = pl.read_parquet(inp / cfg["mapping"]).select(
    "compound_id",
    pl.col("molecule_name").alias("Molecule_Name"),
    pl.col("raw_smiles").alias("SMILES"),
)
assert mapping.height == 750 and mapping["compound_id"].n_unique() == 750
for backend, items in cfg["models"].items():
    assert len(items) == 5
    oofs, tests = [], []
    for item in items:
        summary = json.loads((inp / item["summary"]).read_text())
        assert summary["config"]["backend"] == backend
        assert len(summary["fit_audits"]) == 4 and summary["aux_memberships_verified"]
        oofs.append(pl.read_parquet(inp / item["oof"]))
        tests.append(pl.read_parquet(inp / item["test"]))
    frame = pl.concat(oofs).sort(keys)
    assert frame.select(keys).is_duplicated().sum() == 0
    frames[backend + "_updated"] = frame
    frame.write_parquet(out / f"{backend}-oof.parquet")
    test = pl.concat(tests)
    assert set(test["endpoint"].unique()) == set(base["endpoint"].unique())
    assert test["compound_id"].unique().sort().equals(mapping["compound_id"].sort())
    assert test.height == 750 * 4 * 5 and test.select(keys).is_duplicated().sum() == 0
    counts = test.group_by("compound_id", "endpoint").agg(pl.col("fold").n_unique().alias("n"))
    assert (counts["n"] == 5).all() and test["prediction"].is_finite().all()
    test.write_parquet(out / f"{backend}-test-folds.parquet")
    ensemble = test.group_by("compound_id", "endpoint").agg(pl.col("prediction").mean())
    ensemble.write_parquet(out / f"{backend}-test-ensemble.parquet")
    wide = ensemble.pivot(on="endpoint", index="compound_id", values="prediction")
    regression = mapping.join(wide, on="compound_id", validate="1:1")
    assert (
        regression.height == 750
        and regression.null_count().select(pl.sum_horizontal(pl.all())).item() == 0
    )
    regression.write_parquet(out / f"{backend}-regression.parquet")
rows = []
for method, frame in frames.items():
    for col in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
        assert frame[col].equals(base[col]), (method, col)
    for (ep,), group in frame.group_by("endpoint"):
        for subset in ["all", "neighbors"]:
            selected = group if subset == "all" else group.filter(pl.col("is_neighbor"))
            for fold in [-1, 0, 1, 2, 3, 4]:
                h = selected if fold == -1 else selected.filter(pl.col("fold") == fold)
                y, p, lo, hi = [
                    h[c].to_numpy() for c in ["y_true", "prediction", "conf_low", "conf_high"]
                ]
                assert np.isfinite(p).all()
                rows.append(
                    {
                        "method": method,
                        "endpoint": ep,
                        "subset": subset,
                        "fold": fold,
                        "n": h.height,
                        "ST_RAE": float(
                            official.rae_soft_threshold_absolute_error(
                                y, p, y_true_lower=lo, y_true_upper=hi
                            )
                        ),
                        "MAE": float(np.abs(y - p).mean()),
                        "Spearman": float(spearmanr(y, p).statistic),
                    }
                )
metrics = pl.DataFrame(rows)
metrics.write_csv(out / "metrics.csv")
macro = (
    metrics.filter(pl.col("fold") == -1)
    .group_by("method", "subset")
    .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
    .sort("subset", "method")
)
macro.write_csv(out / "macro.csv")
with psycopg.connect(**cfg["database"]) as conn:
    for r in rows:
        for name in ["ST_RAE", "MAE", "Spearman"]:
            conn.execute(
                "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
                (
                    Path(os.environ["CYP_RUN_DIR"]).name,
                    r["endpoint"],
                    "foundation_replacement:" + r["subset"],
                    str(r["fold"]),
                    r["method"] + ":" + name,
                    r[name],
                ),
            )
write_json(
    out / "summary.json",
    {
        "config": cfg,
        "oof_rows_per_backend": base.height,
        "test_rows_per_backend": 750,
        "models_per_backend": 20,
        "key_label_fold_checks": True,
    },
)
print(macro)
