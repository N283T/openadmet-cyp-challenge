"""Independently verify native trainer cohorts and every outer connectivity mask."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(i):
    return next(inp.glob(f"{i:03d}-*"))


contract = json.loads(read(0).read_text())
bridge = pl.read_parquet(read(1)).filter(pl.col("fold").is_not_null())
records = []
for item in cfg["runs"]:
    table = pl.read_parquet(read(item["table"]))
    assignments = pl.read_parquet(read(item["assignments"]))
    effective = json.loads(read(item["effective"]).read_text())
    heads = effective["heads"]
    model = contract["models"][item["model"]]
    assert heads == model["heads"]
    assert table.height == model["training_rows"] == assignments.height
    assert table["compound_id"].n_unique() == table.height
    joined = table.join(assignments.drop("smiles"), on="compound_id", validate="1:1")
    assert joined.height == table.height
    for h, n in model["label_counts"].items():
        assert np.isfinite(table[h].to_numpy()).sum() == n
    actual_counts = [int(np.isfinite(table[h].to_numpy()).sum()) for h in heads]
    assert actual_counts == effective["global_head_counts"]
    conflicts = (
        joined.group_by("connectivity")
        .agg(pl.col("effective_fold").n_unique().alias("fold_count"))
        .filter(pl.col("fold_count") > 1)
    )
    if conflicts.height:
        joined.join(conflicts, on="connectivity").write_parquet(
            out / f"{item['run']}-connectivity-conflicts.parquet"
        )
    print(item["run"], "connectivity conflicts", conflicts.height, flush=True)
    rows = []
    for fold in range(5):
        train = joined.filter(pl.col("effective_fold") != fold)
        val = joined.filter(pl.col("effective_fold") == fold)
        leaked = set(train["connectivity"]) & set(val["connectivity"])
        if effective["config"]["cv_policy"] == "parent_adapted":
            held = bridge.filter((pl.col("fold") == fold) & ~pl.col("is_anchor"))
            assert not set(train["connectivity"]) & set(held["connectivity"])
            expected = set(table["compound_id"]) & set(held["compound_id"])
            assert expected <= set(val["compound_id"])
            anchors = bridge.filter(pl.col("is_anchor"))
            assert not set(val["connectivity"]) & set(anchors["connectivity"])
        rows.append(
            {
                "fold": fold,
                "cross_fold_connectivity_count": len(leaked),
                "train_rows": train.height,
                "validation_rows": val.height,
                "train_labels": sum(int(np.isfinite(train[h].to_numpy()).sum()) for h in heads),
                "validation_labels": sum(int(np.isfinite(val[h].to_numpy()).sum()) for h in heads),
            }
        )
    for fold, index in item["memberships"].items():
        m = pl.read_parquet(read(index))
        check = m.join(
            assignments.select("compound_id", "effective_fold"), on="compound_id", validate="1:1"
        )
        assert check.height == table.height
        assert (check["training"] == (check["effective_fold"] != int(fold))).all()
        assert (check["validation"] == (check["effective_fold"] == int(fold))).all()
    records.append(
        {
            "run": item["run"],
            "model": item["model"],
            "policy": effective["config"]["cv_policy"],
            "rows": table.height,
            "labels": sum(actual_counts),
            "folds": rows,
            "actual_memberships_verified": list(item["memberships"]),
            "all_five_assignment_masks_verified": True,
            "connectivity_conflicts": conflicts.height,
        }
    )
(out / "audit.json").write_text(json.dumps(records, indent=2))
print(json.dumps(records, indent=2))
