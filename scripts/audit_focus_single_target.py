"""Matched focus endpoint scoring and all-fold exclusion audit."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
new = pl.read_parquet(inp / cfg.get("prediction_file", "000-finetune-oof.parquet"))
base = pl.read_parquet(inp / "001-external-oof.parquet").filter(
    pl.col("endpoint").is_in(new["endpoint"].unique().implode())
)
keys = ["compound_id", "endpoint", "fold"]
new, base = new.sort(keys), base.sort(keys)
for col in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
    assert new[col].equals(base[col]), col
assert new.height == cfg["n"]
assert np.isfinite(new["prediction"].to_numpy()).all()
spec = importlib.util.spec_from_file_location("official", inp / "002-custom_scoring_functions.py")
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
rows = []
for neighbor in [False, True]:
    for fold in [-1, 0, 1, 2, 3, 4]:
        h = new.with_columns(base["prediction"].alias("baseline")).with_columns(
            (0.25 * pl.col("prediction") + 0.75 * pl.col("baseline")).alias("blend25"),
            (0.5 * pl.col("prediction") + 0.5 * pl.col("baseline")).alias("blend50"),
        )
        if neighbor:
            h = h.filter(pl.col("is_neighbor"))
        if fold >= 0:
            h = h.filter(pl.col("fold") == fold)
        if not h.height:
            continue
        for method in ["prediction", "baseline", "blend25", "blend50"]:
            score = scorer.rae_soft_threshold_absolute_error(
                h["y_true"].to_numpy(),
                h[method].to_numpy(),
                y_true_lower=h["conf_low"].to_numpy(),
                y_true_upper=h["conf_high"].to_numpy(),
            )
            rows.append(
                dict(
                    endpoint=cfg["endpoint"],
                    subset="neighbors" if neighbor else "all",
                    fold=fold,
                    method=method,
                    n=h.height,
                    ST_RAE=float(score),
                )
            )
for fold in range(5):
    m = pl.read_parquet(inp / f"{fold + 3:03d}-membership.parquet")
    forbidden = m.filter((pl.col("fold") == fold) & ~pl.col("is_anchor"))
    assert forbidden.filter(pl.col("role") != "outer_holdout").height == 0
    assert set(new.filter(pl.col("fold") == fold)["compound_id"]).issubset(
        set(forbidden["compound_id"])
    )
pl.DataFrame(rows).write_csv(out / "metrics.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(key_label_bound_match=True, outer_membership_exclusion=True, n=new.height), indent=2
    )
)
print(pl.DataFrame(rows).filter((pl.col("fold") == -1) & (pl.col("subset") == "all")))
