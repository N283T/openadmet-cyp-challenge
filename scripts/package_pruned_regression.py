"""Package a target-wise pruned mean from saved full-data regression components."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.tracks import ENDPOINTS, validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoints = ENDPOINTS["regression"]
keys = ["compound_id", "endpoint"]

labels = (
    pl.read_parquet(inp / cfg["measurements"])
    .filter(pl.col("endpoint").is_in(endpoints))
    .select("compound_id", "endpoint", "value")
    .sort(keys)
)
assert labels.height == 6525
reference = pl.read_csv(inp / cfg["reference"])
aliases = (
    pl.read_parquet(inp / cfg["aliases"])
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
test = reference.select("Molecule_Name", "SMILES").join(
    aliases, on="Molecule_Name", how="left", validate="1:1"
)
assert test.height == test["compound_id"].n_unique() == 750
assert test["compound_id"].null_count() == 0
expected = (
    test.select("compound_id").join(pl.DataFrame({"endpoint": endpoints}), how="cross").sort(keys)
)

component_predictions = {}
for component in cfg["components"]:
    prediction = pl.read_parquet(inp / component["predictions"]).sort(keys)
    assert prediction.select(keys).equals(expected)
    assert prediction["split_id"].unique().to_list() == ["full_train"]
    assert prediction["fold"].unique().to_list() == [-1]
    assert prediction["prediction"].is_finite().all()
    membership = pl.read_parquet(inp / component["membership"]).select(labels.columns).sort(keys)
    assert membership.equals(labels), component["name"]
    assert set(membership["compound_id"]).isdisjoint(test["compound_id"])
    component_predictions[component["name"]] = prediction

assert set(cfg["endpoint_components"]) == set(endpoints)
long_parts = []
stats = []
for endpoint in endpoints:
    selected = cfg["endpoint_components"][endpoint]
    assert selected and set(selected).issubset(component_predictions)
    arrays = [
        component_predictions[name].filter(pl.col("endpoint") == endpoint)["prediction"].to_numpy()
        for name in selected
    ]
    prediction = np.mean(np.stack(arrays), axis=0)
    assert np.isfinite(prediction).all()
    part = expected.filter(pl.col("endpoint") == endpoint).with_columns(
        pl.Series("prediction", prediction),
        pl.lit("full_train").alias("split_id"),
        pl.lit(-1).alias("fold"),
    )
    long_parts.append(part)
    stats.append(
        {
            "endpoint": endpoint,
            "component_count": len(selected),
            "components": "+".join(selected),
            "mean": float(prediction.mean()),
            "std": float(prediction.std()),
            "min": float(prediction.min()),
            "max": float(prediction.max()),
        }
    )

long = pl.concat(long_parts).sort(keys).join(test, on="compound_id", validate="m:1")
long.write_parquet(out / "test-predictions.parquet")
wide = (
    long.pivot(on="endpoint", index=["Molecule_Name", "SMILES"], values="prediction")
    .select("SMILES", "Molecule_Name", *endpoints)
    .sort("Molecule_Name")
)
wide.write_parquet(out / "regression.parquet")
wide.write_csv(out / "regression.csv")
validation = validate(out / "regression.parquet", "regression", inp / cfg["reference"])
assert validation["valid"]
assert validate(out / "regression.csv", "regression", inp / cfg["reference"])["valid"]
pl.DataFrame(stats).write_csv(out / "prediction-distributions.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "validation": validation,
            "endpoint_components": cfg["endpoint_components"],
            "selection_rule": cfg["selection_rule"],
            "selection_run": cfg["selection_run"],
            "endpoint_component_count": sum(map(len, cfg["endpoint_components"].values())),
            "previous_endpoint_component_count": 20,
            "calibration": False,
            "deployment": "All-data component predictions; target-wise equal means",
        },
        indent=2,
    )
)
print(json.dumps({"validation": validation, "distributions": stats}, indent=2), flush=True)
