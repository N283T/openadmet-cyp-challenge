"""Verify full-data component memberships and freeze the uncalibrated equal mean."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.tracks import ENDPOINTS, validate

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
labels = (
    pl.read_parquet(inp / cfg["measurements"])
    .filter(pl.col("endpoint").is_in(ENDPOINTS["regression"]))
    .select("compound_id", "endpoint", "value")
    .sort("compound_id", "endpoint")
)
assert labels.height == 6525
ref = pl.read_csv(inp / cfg["reference"])
aliases = (
    pl.read_parquet(inp / cfg["aliases"])
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
test = ref.select("Molecule_Name", "SMILES").join(
    aliases, on="Molecule_Name", how="left", validate="1:1"
)
assert test.height == test["compound_id"].n_unique() == 750
assert test["compound_id"].null_count() == 0
expected = (
    test.select("compound_id")
    .join(pl.DataFrame({"endpoint": ENDPOINTS["regression"]}), how="cross")
    .sort("compound_id", "endpoint")
)
keys = ["compound_id", "endpoint"]
arrays, manifest, stats, pairwise = [], [], [], []
assert len(cfg["components"]) == 5
for c in cfg["components"]:
    assert c["weight"] == 0.2
    pred = pl.read_parquet(inp / c["predictions"]).sort(keys)
    assert pred.select(keys).equals(expected) and pred.height == 3000
    assert pred["split_id"].unique().to_list() == ["full_train"]
    assert pred["fold"].unique().to_list() == [-1]
    if "SMILES" in pred.columns:
        identity = pred.select("compound_id", "Molecule_Name", "SMILES").unique()
        assert identity.sort("compound_id").equals(
            test.select(identity.columns).sort("compound_id")
        )
    membership = pl.read_parquet(inp / c["membership"]).select(labels.columns).sort(keys)
    assert membership.equals(labels), c["name"]
    assert set(membership["compound_id"]).isdisjoint(test["compound_id"])
    values = pred["prediction"].to_numpy()
    assert np.isfinite(values).all()
    arrays.append(values)
    manifest.append(
        dict(
            **c,
            train_labels=membership.height,
            test_predictions=pred.height,
            prediction_sha256=hashlib.sha256((inp / c["predictions"]).read_bytes()).hexdigest(),
            all_original_labels_used=True,
            test_train_overlap=0,
        )
    )
    for ep in ENDPOINTS["regression"]:
        v = pred.filter(pl.col("endpoint") == ep)["prediction"].to_numpy()
        stats.append(
            dict(
                component=c["name"],
                endpoint=ep,
                mean=float(v.mean()),
                std=float(v.std()),
                min=float(v.min()),
                max=float(v.max()),
            )
        )
mean = np.mean(np.stack(arrays), axis=0)
assert np.isfinite(mean).all()
long = (
    expected.with_columns(pl.Series("prediction", mean))
    .join(test, on="compound_id", validate="m:1")
    .with_columns(pl.lit("full_train").alias("split_id"), pl.lit(-1).alias("fold"))
)
long.write_parquet(out / "test-predictions.parquet")
wide = (
    long.pivot(on="endpoint", index=["Molecule_Name", "SMILES"], values="prediction")
    .select("SMILES", "Molecule_Name", *ENDPOINTS["regression"])
    .sort("Molecule_Name")
)
wide.write_parquet(out / "regression.parquet")
wide.write_csv(out / "regression.csv")
validation = validate(out / "regression.parquet", "regression", inp / cfg["reference"])
assert validate(out / "regression.csv", "regression", inp / cfg["reference"])["valid"]
for ep in ENDPOINTS["regression"]:
    mask = expected["endpoint"].to_numpy() == ep
    v = mean[mask]
    stats.append(
        dict(
            component="equal_mean",
            endpoint=ep,
            mean=float(v.mean()),
            std=float(v.std()),
            min=float(v.min()),
            max=float(v.max()),
        )
    )
    for i in range(5):
        for j in range(i + 1, 5):
            pairwise.append(
                dict(
                    endpoint=ep,
                    left=cfg["components"][i]["name"],
                    right=cfg["components"][j]["name"],
                    pearson=float(np.corrcoef(arrays[i][mask], arrays[j][mask])[0, 1]),
                )
            )
pl.DataFrame(stats).write_csv(out / "prediction-distributions.csv")
pl.DataFrame(pairwise).write_csv(out / "component-correlations.csv")
(out / "manifest.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            components=manifest,
            validation=validation,
            calibration=False,
            weights_fitted=False,
            deployment="All-data refit; equal five-component arithmetic mean",
            oof_reference=0.5964141092282933,
            oof_is_development_score_not_test_score=True,
        ),
        indent=2,
    )
)
print(json.dumps(dict(validation=validation, components=manifest, distributions=stats)), flush=True)
