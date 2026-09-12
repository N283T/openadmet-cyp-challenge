"""Compare identical native recipe under parent and published scaffold splits."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import pearsonr, spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


parent, source, reference = [pl.read_parquet(path(i)) for i in range(3)]
reference = reference.filter(
    (pl.col("endpoint") == "CYP2D6_pIC50_direct_inhibition")
    & (pl.col("method") == "caruana_weight_crossfit")
).sort("compound_id")
spec = importlib.util.spec_from_file_location("official", path(3))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
assert len(parent) == 1457 and len(source) == 1493 and len(reference) == 1457
rows = []
parts = []
for name, h in [
    ("native_parent", parent),
    ("native_source_scaffold", source),
    ("previous_global", reference),
]:
    h = h.sort("compound_id")
    assert h["compound_id"].n_unique() == len(h)
    common = h.filter(pl.col("compound_id").is_in(reference["compound_id"].implode()))
    assert common["compound_id"].equals(reference["compound_id"])
    assert np.allclose(common["y_true"], reference["y_true"], rtol=0, atol=1e-9)
    parts.append(common.with_columns(pl.lit(name).alias("method")))
    for cohort, t in [("common_1457", common)] + (
        [("source_all_1493", h)] if len(h) == 1493 else []
    ):
        y, p, lo, hi = [t[c].to_numpy() for c in ["y_true", "prediction", "conf_low", "conf_high"]]
        rows.append(
            dict(
                method=name,
                cohort=cohort,
                n=len(t),
                ST_RAE=float(
                    scorer.rae_soft_threshold_absolute_error(y, p, y_true_lower=lo, y_true_upper=hi)
                ),
                author_local=float(
                    (np.maximum(lo - p, 0) + np.maximum(p - hi, 0)).sum()
                    / np.abs(y - y.mean()).sum()
                ),
                MAE=float(np.abs(y - p).mean()),
                Spearman=float(spearmanr(y, p).statistic),
                Pearson=float(pearsonr(y, p).statistic),
                prediction_mean=float(p.mean()),
                prediction_std=float(p.std()),
            )
        )
pl.concat(parts, how="diagonal_relaxed").write_parquet(out / "comparison-oof.parquet")
pl.DataFrame(rows).write_csv(out / "metrics.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            common_ids=1457,
            source_all=1493,
            same_truth=True,
            folds=(
                "different by design; native source and parent "
                "both select checkpoints on validation"
            ),
        ),
        indent=2,
    )
)
print(pl.DataFrame(rows).write_csv(), flush=True)
