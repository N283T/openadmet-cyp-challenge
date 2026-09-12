"""Score a completed native member against the fixed prior OOF reference."""

import importlib.util
import os
from pathlib import Path

import polars as pl
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


spec = importlib.util.spec_from_file_location("official", path(0))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
ref = pl.read_parquet(path(1)).filter(pl.col("method") == "caruana_weight_crossfit")
member = pl.read_parquet(path(2))
keys = ["compound_id", "endpoint", "fold"]
rows = []
for (endpoint,), t in member.group_by("endpoint"):
    t = t.sort(keys)
    r = ref.filter(pl.col("endpoint") == endpoint).sort(keys)
    assert t.select(keys).equals(r.select(keys))
    for name, h in [("native", t), ("previous_global", r)]:
        y, p, lo, hi = [h[c].to_numpy() for c in ["y_true", "prediction", "conf_low", "conf_high"]]
        rows.append(
            dict(
                method=name,
                endpoint=endpoint,
                n=len(h),
                ST_RAE=float(
                    scorer.rae_soft_threshold_absolute_error(
                        y,
                        p,
                        y_true_lower=lo,
                        y_true_upper=hi,
                    )
                ),
                Spearman=float(spearmanr(y, p).statistic),
            )
        )
result = pl.DataFrame(rows).sort("endpoint", "method")
result.write_csv(out / "metrics.csv")
print(result.write_csv(), flush=True)
