"""Validate and package all ten adopted auxiliary models and keyed predictions."""

import itertools
import json
import os
from pathlib import Path

import polars as pl

from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
selected = json.loads((inp / cfg["selected"]).read_text())
compounds = pl.read_parquet(inp / cfg["compounds"])
assignments = pl.read_parquet(inp / cfg["assignments"])
all_ids = compounds["compound_id"].sort()
models, frames, pairs = [], [], []
for item in cfg["models"]:
    meta = json.loads((inp / item["summary"]).read_text())
    recipe = meta["config"]
    for key in [
        "learning_rate",
        "lr_ratio",
        "encoder_lr_scale",
        "weight_decay",
        "ffn_dropout",
        "batch_norm",
        "scheduler",
        "optimizer",
        "max_epochs",
        "patience",
        "seed",
    ]:
        assert recipe[key] == selected["config"][key], key
    assert recipe["selection_run"] == cfg["selection_run"] and recipe["selection_trial"] == 14
    pair = recipe["excluded_folds"]
    pairs.append(tuple(pair))
    assert recipe["split_run"] == cfg["split_run"]
    members = pl.read_parquet(inp / item["membership"]).join(
        assignments.select("compound_id", "fold"), on="compound_id", validate="1:1"
    )
    assert (members["fold"].is_in(pair) == (members["role"] == "outer_holdout")).all()
    pred = pl.concat([pl.read_parquet(inp / p) for p in item["predictions"]]).sort("compound_id")
    assert pred["compound_id"].equals(all_ids)
    heads = meta["heads"]
    assert (
        len(heads) == 4 and pred.select(pl.all_horizontal(pl.col(heads).is_finite()).all()).item()
    )
    frames.append(
        pred.with_columns(pl.lit(pair[0]).alias("excluded_a"), pl.lit(pair[1]).alias("excluded_b"))
    )
    models.append({"run": item["run"], "excluded_folds": pair})
assert sorted(pairs) == list(itertools.combinations(range(5), 2))
joined = pl.concat(frames)
joined.write_parquet(out / "nested-predictions.parquet")
oof = []
for fold in range(5):
    fold_predictions = joined.filter(
        (pl.col("excluded_a") == fold) | (pl.col("excluded_b") == fold)
    )
    p = fold_predictions.group_by("compound_id").agg(
        pl.col(heads).mean(), pl.len().alias("n_models")
    )
    assert (p["n_models"] == 4).all()
    oof.append(
        assignments.filter(pl.col("fold") == fold)
        .select("compound_id", "fold")
        .join(p, on="compound_id", validate="1:1")
    )
pl.concat(oof).with_columns(pl.lit(cfg["split_run"]).alias("split_id")).write_parquet(
    out / "oof-log2fc.parquet"
)
test = compounds.select("compound_id").join(
    assignments.select("compound_id"), on="compound_id", how="anti"
)
assert test.height == 750
test.join(
    joined.group_by("compound_id").agg(pl.col(heads).mean()), on="compound_id", validate="1:1"
).write_parquet(out / "test-log2fc.parquet")
write_json(
    out / "registry.json",
    {
        "track": "regression",
        "selection_run": cfg["selection_run"],
        "selection_trial": 14,
        "split_run": cfg["split_run"],
        "models": models,
        "package_run": Path(os.environ["CYP_RUN_DIR"]).name,
        "adoption": "User requested replacement with the Rustuna-selected recipe",
        "evaluation_scope": "Fold 0 excluded from HPO; other folds evaluated after selection",
        "recipe": {
            k: v
            for k, v in selected["config"].items()
            if k not in ["outer_fold", "excluded_folds", "tuning_only"]
        },
    },
)
write_json(
    out / "summary.json",
    {
        "models": 10,
        "compounds_per_model": compounds.height,
        "oof_compounds": assignments.height,
        "test_compounds": test.height,
        "pair_memberships_verified": True,
    },
)
