"""Independently reconstruct and score saved global Caruana predictions."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def input_file(index):
    matches = sorted(inp.glob(f"{index:03d}-*"))
    if not matches:
        matches = sorted(inp.glob(f"{index}-*"))
    assert len(matches) == 1, (index, [p.name for p in matches])
    return matches[0]


pool = pl.read_parquet(input_file(0))
predictions = pl.read_parquet(input_file(1))
weights = json.loads(input_file(2).read_text())
spec = importlib.util.spec_from_file_location("official", input_file(3))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
keys = ["compound_id", "endpoint", "fold"]
canonical = keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]
pool = pool.sort(keys)
assert pool.height == 6384
assert not any(pool[c].null_count() for c in canonical)
assert not pool.select(keys).is_duplicated().any()
assert set(pool["fold"]) == set(range(5))
assert all(pool[c].is_finite().all() for c in ["y_true", "conf_low", "conf_high"])
assert (pool["conf_low"] <= pool["conf_high"]).all()
assert set(pool["endpoint"]) == {
    f"CYP{ep}_pIC50_direct_inhibition" for ep in ["1A2", "2C9", "2D6", "3A4"]
}
expected_n = {"1A2": 1372, "2C9": 1255, "2D6": 1457, "3A4": 2300}
for short, n in expected_n.items():
    assert pool.filter(pl.col("endpoint") == f"CYP{short}_pIC50_direct_inhibition").height == n
methods = sorted(predictions["method"].unique())
assert {"caruana_full_fit", "caruana_weight_crossfit", "latest"}.issubset(methods)
for method in methods:
    table = predictions.filter(pl.col("method") == method).sort(keys)
    assert table.height == 6384 and not table.select(keys).is_duplicated().any(), method
    for column in canonical:
        assert table[column].equals(pool[column]), (method, column)
    assert table["prediction"].null_count() == 0, method
    assert table["prediction"].is_finite().all(), method

weight_keys = set()
reconstructions = []
for record in weights:
    endpoint, method, held = record["endpoint"], record["method"], record["held_fold"]
    assert method in {"caruana_full_fit", "caruana_weight_crossfit"}
    assert (held is None) == (method == "caruana_full_fit")
    assert held is None or held in range(5)
    if "selection_fold_ids" in record:
        expected_folds = list(range(5)) if held is None else [f for f in range(5) if f != held]
        assert record["selection_fold_ids"] == expected_folds, record
    key = endpoint, method, held
    assert key not in weight_keys, key
    weight_keys.add(key)
    selected = record["weights"]
    assert selected and all(np.isfinite(v) and v >= 0 for v in selected.values()), key
    assert np.isclose(sum(selected.values()), 1, atol=1e-12, rtol=0), key
    mask = pl.col("endpoint") == endpoint
    if held is not None:
        mask &= pl.col("fold") == held
    source = pool.filter(mask).sort(keys)
    target = predictions.filter(mask & (pl.col("method") == method)).sort(keys)
    assert source.height > 0 and source.select(keys).equals(target.select(keys))
    rebuilt = np.zeros(source.height)
    for column, weight in selected.items():
        assert column.startswith("m") and column[1:].isdigit() and column in source.columns
        values = source[column].to_numpy()
        assert np.isfinite(values).all(), (key, column)
        rebuilt += weight * values
    difference = float(np.max(np.abs(rebuilt - target["prediction"].to_numpy())))
    assert difference < 1e-10, (key, difference)
    reconstructions.append(
        {
            "endpoint": endpoint,
            "method": method,
            "held_fold": held,
            "n": source.height,
            "max_absolute_error": difference,
            "selected_candidates": len(selected),
            "weight_sum": sum(selected.values()),
        }
    )
expected = {(ep, "caruana_full_fit", None) for ep in pool["endpoint"].unique()}
expected |= {
    (ep, "caruana_weight_crossfit", f) for ep in pool["endpoint"].unique() for f in range(5)
}
assert weight_keys == expected, (weight_keys - expected, expected - weight_keys)

rows = []
for (method, endpoint), table in predictions.group_by("method", "endpoint"):
    for subset in ["all", "neighbors"]:
        sub = table if subset == "all" else table.filter(pl.col("is_neighbor"))
        for fold in [-1, 0, 1, 2, 3, 4]:
            group = sub if fold == -1 else sub.filter(pl.col("fold") == fold)
            if not group.height:
                continue
            y, pred = (group[c].to_numpy() for c in ["y_true", "prediction"])
            score = float(
                scorer.rae_soft_threshold_absolute_error(
                    y,
                    pred,
                    y_true_lower=group["conf_low"].to_numpy(),
                    y_true_upper=group["conf_high"].to_numpy(),
                )
            )
            rows.append(
                dict(
                    method=method,
                    endpoint=endpoint,
                    subset=subset,
                    fold=fold,
                    n=group.height,
                    ST_RAE=score,
                    MAE=float(np.abs(y - pred).mean()),
                    Spearman=float(spearmanr(y, pred).statistic),
                )
            )
metrics = pl.DataFrame(rows).sort("method", "subset", "endpoint", "fold")
metrics.write_csv(out / "metrics.csv")
macro = (
    metrics.filter(pl.col("fold") == -1)
    .group_by("method", "subset")
    .agg(pl.col("ST_RAE", "MAE", "Spearman").mean(), pl.len().alias("endpoints"))
    .sort("subset", "ST_RAE")
)
assert (macro["endpoints"] == 4).all()
macro.write_csv(out / "macro.csv")
latest = macro.filter((pl.col("method") == "latest") & (pl.col("subset") == "all"))["ST_RAE"].item()
assert np.isclose(latest, 0.5869306288697727, atol=1e-10, rtol=0), latest
(out / "source-weights.json").write_text(json.dumps(weights, indent=2))
(out / "verification.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            matched_rows=6384,
            methods=methods,
            all_predictions_finite=True,
            canonical_labels_exact=True,
            weights_verified=True,
            reconstructions=reconstructions,
            latest_ST_RAE=latest,
            qualification=(
                "Full-fit score reuses labels for weight selection. Weight-crossfit excludes "
                "held-fold labels from weight fitting only; underlying candidate training "
                "is not newly nested."
            ),
        ),
        indent=2,
    )
)
print(macro)
