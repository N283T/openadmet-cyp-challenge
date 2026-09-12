"""Matched paired readout and fixed ensemble-addition assessment."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
keys = ["compound_id", "endpoint", "fold"]
frames = [pl.read_parquet(next(inp.glob(f"{i:03d}-*"))).sort(keys) for i in range(3)]
b = frames[0]
for f in frames:
    assert f.height == 6384 and f.select(keys + ["y_true", "conf_low", "conf_high"]).equals(
        b.select(keys + ["y_true", "conf_low", "conf_high"])
    )
b = b.with_columns(
    frames[1]["prediction"].alias("old"),
    frames[2]["prediction"].alias("new"),
    pl.col("prediction").alias("final"),
).with_columns(
    ((5 * pl.col("final") + pl.col("new")) / 6).alias("equal_six"),
    ((pl.col("final") + pl.col("new")) / 2).alias("equal_pair"),
)
rows = []
for (ep,), g in b.group_by("endpoint"):
    y, lo, hi = (g[c].to_numpy() for c in ["y_true", "conf_low", "conf_high"])
    den = np.mean(np.maximum(lo - y.mean(), 0) + np.maximum(y.mean() - hi, 0))
    for m in ["old", "new", "final", "equal_six", "equal_pair"]:
        p = g[m].to_numpy()
        rows.append(
            dict(
                endpoint=ep,
                method=m,
                n=g.height,
                st_rae=float(np.mean(np.maximum(lo - p, 0) + np.maximum(p - hi, 0)) / den),
                mae=float(np.abs(p - y).mean()),
                spearman=float(spearmanr(y, p).statistic),
                residual_corr_final=float(np.corrcoef(p - y, g["final"].to_numpy() - y)[0, 1]),
            )
        )
pl.DataFrame(rows).write_csv(out / "comparison.csv")
b.write_parquet(out / "aligned-oof.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "rows": b.height,
            "fit_weights": False,
            "equal_six": "(5 * submitted final + new) / 6",
            "equal_pair": "(submitted final + new) / 2",
        },
        indent=2,
    )
)
print(pl.DataFrame(rows), flush=True)
