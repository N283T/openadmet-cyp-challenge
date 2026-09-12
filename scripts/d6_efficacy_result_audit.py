"""Keyed paired efficacy audit and prespecified ensemble additions; no weight search."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import psycopg
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(index):
    return next(inp.glob(f"{index:03d}-*"))


new, latest, submitted, memberships, assignments = [pl.read_parquet(path(i)) for i in range(5)]
latest = latest.filter(pl.col("endpoint") == "CYP2D6_pIC50_direct_inhibition")
submitted = submitted.filter(pl.col("endpoint") == "CYP2D6_pIC50_direct_inhibition")
spec = importlib.util.spec_from_file_location("official", path(5))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
keys = ["compound_id", "endpoint", "fold"]
latest, submitted = latest.sort(keys), submitted.sort(keys)
assert latest.height == submitted.height == 1457
for col in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
    assert latest[col].equals(submitted[col]), col
for (_ep, fold), m in memberships.group_by("endpoint", "outer_fold"):
    m = m.join(assignments.select("compound_id", "fold"), on="compound_id", validate="1:1")
    expected = (m["fold"].to_numpy() == fold) & ~m["is_anchor"].to_numpy()
    assert np.array_equal(m["validation"].to_numpy(), expected)
predictions = {
    "latest": latest["prediction"].to_numpy(),
    "submitted": submitted["prediction"].to_numpy(),
}
for method in sorted(new["method"].unique()):
    h = new.filter(pl.col("method") == method).sort(keys)
    assert h.height == 1457 and not h.select(keys).is_duplicated().any()
    for col in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
        assert latest[col].equals(h[col]), col
    p = h["prediction"].to_numpy()
    assert np.isfinite(p).all()
    predictions[method] = p
    for base in ["latest", "submitted"]:
        for weight in [0.25, 0.5]:
            predictions[f"{base}__{method}__{int(100 * weight)}pct"] = (1 - weight) * predictions[
                base
            ] + weight * p
aligned = latest.select(keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]).with_columns(
    [pl.Series(name, values) for name, values in predictions.items()]
)
rows = []
for (ep,), e in aligned.group_by("endpoint"):
    for subset in ["all", "neighbors"]:
        g = e if subset == "all" else e.filter(pl.col("is_neighbor"))
        for fold in [-1, 0, 1, 2, 3, 4]:
            h = g if fold == -1 else g.filter(pl.col("fold") == fold)
            y, lo, hi = (h[c].to_numpy() for c in ["y_true", "conf_low", "conf_high"])
            for method in predictions:
                p = h[method].to_numpy()
                rows.append(
                    dict(
                        endpoint=ep,
                        subset=subset,
                        fold=fold,
                        n=h.height,
                        method=method,
                        ST_RAE=float(
                            scorer.rae_soft_threshold_absolute_error(
                                y, p, y_true_lower=lo, y_true_upper=hi
                            )
                        ),
                        MAE=float(np.abs(y - p).mean()),
                        Spearman=float(spearmanr(y, p).statistic),
                    )
                )
metrics = pl.DataFrame(rows)
macro = (
    metrics.filter(pl.col("fold") == -1)
    .group_by("method", "subset")
    .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
    .sort("subset", "ST_RAE")
)
metrics.write_csv(out / "metrics.csv")
macro.write_csv(out / "macro.csv")
aligned.write_parquet(out / "aligned-oof.parquet")
with psycopg.connect(**cfg["database"]) as conn:
    for row in rows:
        conn.execute(
            "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
            (
                Path(os.environ["CYP_RUN_DIR"]).name,
                row["endpoint"],
                "efficacy:" + row["subset"],
                str(row["fold"]),
                row["method"] + ":ST_RAE",
                row["ST_RAE"],
            ),
        )
(out / "verification.json").write_text(
    json.dumps(
        dict(
            matched_per_method=1457,
            exact_keys_labels_bounds=True,
            membership_verified=True,
            fixed_weights=[0.25, 0.5],
            weight_search=False,
            interpretation="Development comparison, no selection-adjusted independent estimate",
        ),
        indent=2,
    )
)
print(macro)
