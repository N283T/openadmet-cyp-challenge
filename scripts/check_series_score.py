"""Decompose series ST-RAE and independently verify keyed predictions."""

import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
oof = pl.read_parquet(next(inp.glob("*-oof.parquet")))
labels = pl.read_parquet(next(inp.glob("*-measurements.parquet")))
seeds = set(pl.read_csv(next(inp.glob("*-hit-seeds.csv")))["compound_id"])
rows = []
for (strategy, ep, model), f in oof.group_by("strategy", "endpoint", "model"):
    ref = labels.filter(pl.col("endpoint") == ep).select(
        "compound_id", pl.col("value").alias("label")
    )
    f = f.join(ref, on="compound_id", validate="1:1", maintain_order="left")
    assert np.allclose(f["label"], f["y_true"])
    for subset in ["all", "neighbors"]:
        g = (
            f
            if subset == "all"
            else f.filter((pl.col("hit_origins") != "") & ~pl.col("compound_id").is_in(list(seeds)))
        )
        y, p, lo, hi = [g[c].to_numpy() for c in ["y_true", "prediction", "conf_low", "conf_high"]]
        err = np.maximum(lo - p, 0) + np.maximum(p - hi, 0)
        den = np.maximum(lo - y.mean(), 0) + np.maximum(y.mean() - hi, 0)
        rows.append(
            dict(
                strategy=strategy,
                endpoint=ep,
                model=model,
                subset=subset,
                n=len(y),
                y_mean=float(y.mean()),
                y_sd=float(y.std()),
                prediction_mean=float(p.mean()),
                bias=float((p - y).mean()),
                MAE=float(np.abs(p - y).mean()),
                soft_MAE=float(err.mean()),
                denominator=float(den.mean()),
                ST_RAE=float(err.sum() / den.sum()),
                ci_width=float((hi - lo).mean()),
            )
        )
pl.DataFrame(rows).sort("strategy", "subset", "endpoint", "model").write_csv(
    out / "decomposition.csv"
)
write_json(
    out / "summary.json",
    dict(
        key_alignment_verified=True,
        definition="ST_RAE = mean interval error / mean interval error of subset label mean",
    ),
)
print(pl.DataFrame(rows).filter((pl.col("strategy") == "series") & (pl.col("model") == "tabpfn")))
