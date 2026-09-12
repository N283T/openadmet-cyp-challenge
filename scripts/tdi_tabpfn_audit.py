"""Audit complete keyed probabilities and summarize threshold-independent metrics."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from tdi_tabpfn_frozen import choose_threshold

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
oof = pl.read_parquet(inp / cfg["oof"])
inner = pl.read_parquet(inp / cfg["inner"])
test = pl.read_parquet(inp / cfg["test"])
fold_test = pl.read_parquet(inp / cfg["fold_test"])
labels = pl.read_parquet(inp / cfg["labels"])
membership_counts = {"outer": 0, "inner": 0}
for member in cfg["memberships"]:
    table = pl.read_parquet(inp / member["file"])
    fold = member["fold"]
    outer_mask = table["fold"] == fold
    if member["kind"] == "outer":
        assert table["validation"].equals(outer_mask.rename("validation"))
        assert table["training"].equals((~outer_mask).rename("training"))
    else:
        inner_mask = table["fold"] == (fold + 1) % 5
        assert table["inner_validation"].equals(inner_mask.rename("inner_validation"))
        assert table["training"].equals((~(outer_mask | inner_mask)).rename("training"))
    membership_counts[member["kind"]] += 1
assert membership_counts == {"outer": 30, "inner": 30}
validation = pl.read_parquet(inp / cfg["inner_validation"])
assert validation.height == 15243
for (method, endpoint, fold), table in validation.group_by("method", "endpoint", "outer_fold"):
    assert set(table["inner_fold"]) == {(fold + 1) % 5}
    expected = labels.filter((pl.col("endpoint") == endpoint) & (pl.col("fold") == (fold + 1) % 5))
    assert set(table["compound_id"]) == set(expected["compound_id"])
    threshold = choose_threshold(table["y_true"].to_numpy(), table["probability"].to_numpy())
    scored = inner.filter(
        (pl.col("method") == method) & (pl.col("endpoint") == endpoint) & (pl.col("fold") == fold)
    )
    np.testing.assert_array_equal(scored["inner_threshold"], np.full(scored.height, threshold))
    np.testing.assert_array_equal(
        scored["prediction_inner"], (scored["probability"].to_numpy() >= threshold).astype(int)
    )
keys = ["compound_id", "endpoint", "method", "fold"]
assert oof.height == inner.height == 15243
assert oof.unique(keys).height == oof.height
assert test.height == 4500 and test.unique(["compound_id", "endpoint", "method"]).height == 4500
assert fold_test.height == 22500
paired = oof.join(
    inner.select(*keys, pl.col("probability").alias("inner_probability"), "prediction_inner"),
    on=keys,
    validate="1:1",
)
assert paired.height == oof.height
np.testing.assert_array_equal(paired["probability"], paired["inner_probability"])
rows = []
for (method, endpoint), table in paired.group_by("method", "endpoint"):
    expected = labels.filter(pl.col("endpoint") == endpoint)
    assert table.height == (1497 if endpoint == "CYP2D6_is_TDI" else 3584)
    assert set(table.select("compound_id", "fold").iter_rows()) == set(
        expected.select("compound_id", "fold").iter_rows()
    )
    comparison = table.join(
        expected.select("compound_id", pl.col("value").alias("label")),
        on="compound_id",
        validate="1:1",
    )
    np.testing.assert_array_equal(comparison["y_true"], comparison["label"])
    assert table["probability"].is_between(0, 1).all()
    y, p = table["y_true"].to_numpy(), table["probability"].to_numpy()
    rows.append(
        dict(
            method=method,
            endpoint=endpoint,
            n=len(y),
            positive_rate=float(y.mean()),
            average_precision=average_precision_score(y, p),
            roc_auc=roc_auc_score(y, p),
            mcc_05=matthews_corrcoef(y, p >= 0.5),
            mcc_inner=matthews_corrcoef(y, table["prediction_inner"]),
            probability_mean=float(p.mean()),
            probability_min=float(p.min()),
            probability_max=float(p.max()),
            positive_predictions_05=int((p >= 0.5).sum()),
        )
    )
rebuilt = fold_test.group_by("compound_id", "endpoint", "method").agg(
    pl.col("probability").mean().alias("rebuilt"), pl.len().alias("count")
)
paired_test = test.join(rebuilt, on=["compound_id", "endpoint", "method"], validate="1:1")
assert (paired_test["count"] == 5).all()
np.testing.assert_allclose(paired_test["probability"], paired_test["rebuilt"], atol=1e-12, rtol=0)
pl.DataFrame(rows).sort("endpoint", "method").write_csv(out / "metrics.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            oof_rows=15243,
            test_rows=4500,
            probabilities_unchanged=True,
            canonical_keys_labels_exact=True,
            test_mean_max_error=float(
                (paired_test["probability"] - paired_test["rebuilt"]).abs().max()
            ),
        ),
        indent=2,
    )
)
