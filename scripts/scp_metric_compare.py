"""Compare exact source metric implementations on identical saved OOF predictions."""

import ast
import importlib.util
import json
import logging
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


predictions = pl.read_parquet(path(0))
spec = importlib.util.spec_from_file_location("saved_official", path(1))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
# Execute only the two unmodified pure metric functions from the pinned source;
# this avoids importing Workbench's AWS/application dependencies.
tree = ast.parse(path(2).read_text())
selected = [
    node
    for node in tree.body
    if isinstance(node, ast.FunctionDef)
    and node.name in ["soft_threshold_error", "soft_threshold_rae"]
]
assert len(selected) == 2
namespace = {"np": np, "log": logging.getLogger("metric-comparison")}
exec(compile(ast.Module(body=selected, type_ignores=[]), str(path(2)), "exec"), namespace)
rows = []
for (method, endpoint), h in predictions.group_by("method", "endpoint"):
    y, p, lo, hi = [h[c].to_numpy() for c in ["y_true", "prediction", "conf_low", "conf_high"]]
    saved = float(scorer.rae_soft_threshold_absolute_error(y, p, y_true_lower=lo, y_true_upper=hi))
    source_default = namespace["soft_threshold_rae"](y, p, lo, hi)
    source_soft = namespace["soft_threshold_rae"](y, p, lo, hi, soft_baseline=True)
    assert np.isclose(saved, source_soft, rtol=0, atol=1e-12)
    numerator = float((np.maximum(lo - p, 0) + np.maximum(p - hi, 0)).sum())
    plain_denominator = float(np.abs(y - y.mean()).sum())
    soft_denominator = float((np.maximum(lo - y.mean(), 0) + np.maximum(y.mean() - hi, 0)).sum())
    assert np.isclose(source_default, saved * soft_denominator / plain_denominator)
    rows.append(
        dict(
            method=method,
            endpoint=endpoint,
            n=len(h),
            saved_scorer=saved,
            author_default=source_default,
            author_soft_baseline=source_soft,
            numerator=numerator,
            plain_denominator=plain_denominator,
            soft_denominator=soft_denominator,
            ratio=soft_denominator / plain_denominator,
        )
    )
r = pl.DataFrame(rows).sort("method", "endpoint")
r.write_csv(out / "metric-comparison.csv")
r.group_by("method").agg(pl.col("saved_scorer").mean(), pl.col("author_default").mean()).write_csv(
    out / "macro-comparison.csv"
)
(out / "verification.json").write_text(
    json.dumps(
        dict(
            same_predictions=True,
            unchanged_source_functions=True,
            saved_equals_source_soft=True,
            live_leaderboard_formula="not established by this run",
        ),
        indent=2,
    )
)
print(r.filter(pl.col("method") == "caruana_weight_crossfit").write_csv(), flush=True)
