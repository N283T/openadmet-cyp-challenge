"""Independently verify deployment membership and frozen selection transfer."""

import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
source = json.loads((inp / cfg["source_selection"]).read_text())
actual = json.loads((inp / cfg["actual_selection"]).read_text())
labels = (
    pl.read_parquet(inp / cfg["measurements"])
    .filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))
    .select("compound_id", "endpoint", "value")
    .sort("compound_id", "endpoint")
)
member = (
    pl.read_parquet(inp / cfg["membership"])
    .select("compound_id", "endpoint", "value")
    .sort("compound_id", "endpoint")
)
assert labels.equals(member)
pred = pl.read_parquet(inp / cfg["predictions"])
test = pl.read_parquet(inp / cfg["test_mapping"]).select("compound_id").unique()
assert labels.height == 6525 and test.height == 750 and pred.height == 3000
assert pred.select("compound_id", "endpoint").n_unique() == 3000
assert np.isfinite(pred["prediction"].to_numpy()).all()
assert member.join(test, on="compound_id", how="inner").is_empty()
checks = []
for ep in sorted(labels["endpoint"].unique()):
    expected = defaultdict(float)
    for fold in source["selections"]:
        if fold["endpoint"] == ep:
            for c in fold["selected"]:
                expected[c["name"]] += c["weight"] / 5
    chosen = next(x for x in actual["endpoints"] if x["endpoint"] == ep)
    weights = {c["name"]: c["weight"] for c in chosen["selected"]}
    assert weights == dict(expected)
    assert np.isclose(sum(weights.values()), 1)
    n_train = labels.filter(pl.col("endpoint") == ep).height
    assert chosen["n_train"] == n_train
    assert all(c["n_train"] == n_train for c in chosen["selected"])
    ep_test = pred.filter(pl.col("endpoint") == ep)
    assert ep_test.select("compound_id").sort("compound_id").equals(test.sort("compound_id"))
    checks.append(
        dict(
            endpoint=ep,
            n_train=n_train,
            n_test=ep_test.height,
            n_selected=len(weights),
            sum_weights=sum(weights.values()),
        )
    )
(out / "audit.json").write_text(
    json.dumps(
        dict(
            passed=True,
            checks=checks,
            all_endpoint_labels_retained=True,
            test_excluded_from_training=True,
            weights_match_saved_inner_selections=True,
        ),
        indent=2,
    )
)
print(json.dumps(checks), flush=True)
