"""Assemble fold-safe auxiliary predictions and audit all model training memberships."""

import itertools
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
reg = json.loads((inp / cfg["registry"]).read_text())
readout = json.loads((inp / cfg["readout_registry"]).read_text())
a = pl.read_parquet(inp / cfg["parent_assignments"]).sort("compound_id")
original = pl.read_parquet(inp / cfg["aux_assignments"])
print("parent_columns", a.columns, "aux_columns", original.columns, flush=True)
assert a.height == a["compound_id"].n_unique() == 6146
assert set(a["fold"]) == set(range(5))
assert readout["auxiliary_package_run"] == reg["package_run"]
if "strategy" in original.columns:
    original = original.filter(pl.col("strategy") == reg["recipe"]["outer_strategy"])
check = a.select("compound_id", "fold").join(
    original.select("compound_id", pl.col("fold").alias("aux_fold")),
    on="compound_id",
    validate="1:1",
)
assert check.height == a.height and (check["fold"] == check["aux_fold"]).all()
assert a.group_by("group").agg(pl.col("fold").n_unique().alias("n"))["n"].max() == 1
assert a.group_by("connectivity").agg(pl.col("fold").n_unique().alias("n"))["n"].max() == 1
models, audits = {}, []
head_order = None
for item in cfg["models"]:
    meta = json.loads((inp / item["summary"]).read_text())
    pair = tuple(sorted(meta["config"]["excluded_folds"]))
    assert list(pair) == item["excluded_folds"]
    assert meta["config"]["split_run"] == reg["split_run"]
    heads = meta["heads"]
    head_order = heads if head_order is None else head_order
    assert heads == head_order and len(heads) == 4
    members = pl.read_parquet(inp / item["membership"]).join(
        a.select("compound_id", "fold", "is_anchor"), on="compound_id", validate="1:1"
    )
    assert members.height == sum(meta["counts"].values())
    held = members["fold"].is_in(pair)
    assert (held == (members["role"] == "outer_holdout")).all()
    assert set(members["role"]) <= {
        "outer_holdout",
        "inner_train",
        "inner_validation",
        "train",
        "validation",
    }
    pred = pl.concat(pl.read_parquet(inp / p) for p in item["predictions"])
    assert pred.height == pred["compound_id"].n_unique() == 6896
    v = (
        a.select("compound_id")
        .join(pred, on="compound_id", validate="1:1", maintain_order="left")
        .select(heads)
        .to_numpy()
        .astype(np.float32)
    )
    assert np.isfinite(v).all()
    assert pair not in models
    models[pair] = v
    audits.append(
        {
            "run": item["run"],
            "excluded_folds": list(pair),
            "members": members.height,
            "heldout": int(held.sum()),
            "heldout_anchors": members.filter(held & members["is_anchor"]).height,
            "heads": heads,
            "prediction_rows": pred.height,
        }
    )
assert set(models) == set(itertools.combinations(range(5), 2))
frames = []
ff = a["fold"].to_numpy()
for outer in range(5):
    # This follows the former ensemble_components.py float32 accumulation exactly.
    aux = np.zeros((a.height, 4), dtype=np.float32)
    for pair, v in models.items():
        if outer not in pair:
            continue
        other = next(f for f in pair if f != outer)
        aux[ff == other] = v[ff == other]
        aux[ff == outer] += v[ff == outer] / 4
    # Independent recipient-wise reconstruction tests every matrix entry, including anchors.
    reference = np.zeros_like(aux)
    for recipient in range(5):
        mask = ff == recipient
        if recipient != outer:
            reference[mask] = models[tuple(sorted((outer, recipient)))][mask]
        else:
            stack = np.stack([v[mask] for pair, v in models.items() if outer in pair])
            reference[mask] = stack.sum(axis=0, dtype=np.float32) / 4
    assert np.array_equal(aux, reference)
    f = a.with_columns(
        pl.lit(outer).alias("outer_fold"),
        pl.lit(readout["split_run"]).alias("split_id"),
        ((pl.col("fold") == outer) & ~pl.col("is_anchor")).alias("is_outer_validation"),
        pl.when(pl.col("fold") == outer).then(4).otherwise(1).alias("aux_model_count"),
        *[pl.Series(f"aux_{i}", aux[:, i]) for i in range(4)],
    )
    frames.append(f)
result = pl.concat(frames)
assert result.height == result.select("compound_id", "outer_fold").n_unique() == 30730
result.write_parquet(out / "nested-log2fc.parquet")
(out / "membership-audit.json").write_text(json.dumps(audits, indent=2))
summary = dict(
    rows=result.height,
    compounds=a.height,
    outer_folds=5,
    heads={f"aux_{i}": h for i, h in enumerate(head_order)},
    split_run=readout["split_run"],
    aux_split_run=reg["split_run"],
    registry_package_run=reg["package_run"],
    matrix_check="Bitwise equal to independent recipient-wise reconstruction in all five folds",
    prior_check="Same routing as ensemble_components.py; no saved matrix for direct comparison",
    train_features_exclude_outer_and_recipient=True,
    anchors_also_excluded_in_auxiliary=True,
)
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary), flush=True)
