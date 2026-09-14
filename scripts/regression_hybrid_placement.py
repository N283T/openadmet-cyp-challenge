"""Build validated hybrid placement scenarios from the current simplified vector."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.tracks import ENDPOINTS, validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    path = inp / cfg[name]
    return pl.read_csv(path) if path.suffix == ".csv" else pl.read_parquet(path)


raw = read("predictions").select("compound_id", "endpoint", "prediction").sort(
    "compound_id", "endpoint"
)
aliases = (
    read("aliases")
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
reference_path = inp / cfg["reference"]
reference = pl.read_csv(reference_path)
identity = reference.select("Molecule_Name", "SMILES").join(
    aliases, on="Molecule_Name", validate="1:1"
)
assert raw.height == 3000 and raw.select("compound_id", "endpoint").n_unique() == 3000
assert identity.height == identity["compound_id"].n_unique() == 750


def place(values: np.ndarray, mean: float, sd: float) -> np.ndarray:
    return mean + sd * (values - values.mean()) / values.std(ddof=1)


scenario_rows = []
distribution_rows = []
manifest = []
for scenario in cfg["scenarios"]:
    frames = []
    for (endpoint,), group in raw.group_by("endpoint", maintain_order=True):
        values = group["prediction"].to_numpy().astype(float)
        target = cfg["targets"].get(endpoint)
        strength = scenario["strengths"].get(endpoint, 0.0)
        if target is None or strength == 0:
            prediction = values.copy()
        else:
            fully_placed = place(values, target["mean"], target["sd"])
            prediction = values + strength * (fully_placed - values)
        frames.append(group.with_columns(pl.Series("prediction", prediction)))
        distribution_rows.append(
            {
                "scenario": scenario["name"],
                "endpoint": endpoint,
                "strength": strength,
                "raw_mean": float(values.mean()),
                "raw_sd": float(values.std(ddof=1)),
                "new_mean": float(prediction.mean()),
                "new_sd": float(prediction.std(ddof=1)),
                "mean_abs_change": float(np.abs(prediction - values).mean()),
                "minimum": float(prediction.min()),
                "maximum": float(prediction.max()),
            }
        )
    long = pl.concat(frames).sort("compound_id", "endpoint")
    long = long.with_columns(pl.lit(scenario["name"]).alias("scenario"))
    scenario_rows.append(long)
    wide = (
        long.join(identity, on="compound_id", validate="m:1")
        .pivot(on="endpoint", index=["Molecule_Name", "SMILES"], values="prediction")
        .select("SMILES", "Molecule_Name", *ENDPOINTS["regression"])
        .sort("Molecule_Name")
    )
    destination = out / f"{scenario['name']}.csv"
    wide.write_csv(destination)
    checked = validate(destination, "regression", reference_path)
    assert checked["valid"], checked
    manifest.append(
        {
            "scenario": scenario["name"],
            "file": destination.name,
            "sha256": hashlib.sha256(destination.read_bytes()).hexdigest(),
            "validation": checked,
        }
    )

pl.concat(scenario_rows).write_parquet(out / "scenario-predictions.parquet")
pl.DataFrame(distribution_rows).sort("scenario", "endpoint").write_csv(
    out / "distributions.csv"
)
(out / "manifest.json").write_text(
    json.dumps(
        {
            "source": "Current simplified regression test predictions",
            "formula": "raw + strength * (full_affine_placement - raw)",
            "targets": cfg["targets"],
            "files": manifest,
            "score_claim": "None; source-assumption scenarios only",
        },
        indent=2,
    )
    + "\n"
)
print(pl.DataFrame(distribution_rows).sort("scenario", "endpoint"), flush=True)
