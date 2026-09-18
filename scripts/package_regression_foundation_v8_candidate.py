"""Package the v8 target-specific flat banks as a regression candidate."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

from cyp_challenge.tracks import validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

aliases = (
    pl.read_parquet(inp / cfg["aliases"])
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
    .sort("compound_id")
)
reference = pl.read_csv(inp / cfg["reference"])
test = aliases.join(reference, on="Molecule_Name", validate="1:1", maintain_order="left")
if test.height != test["compound_id"].n_unique() or test.height != 750:
    raise ValueError("Invalid blinded-test identity table")

atomic = pl.read_parquet(inp / cfg["prior_atomic"]).select(
    "compound_id", "endpoint", "member", "prediction"
)
new = pl.read_parquet(inp / cfg["new_atomic"]).select(
    "compound_id", "endpoint", pl.col("method").alias("source_method"), "prediction"
)
rename_rows = []
for source, member in cfg["new_method_names"].items():
    rows = new.filter(pl.col("source_method") == source).select(
        "compound_id",
        "endpoint",
        pl.lit(member).alias("member"),
        "prediction",
    )
    if rows.height != 750:
        raise ValueError(f"Missing new full-fit prediction: {source}")
    rename_rows.append(rows)

for source in cfg["supplemental_predictions"]:
    frame = pl.read_parquet(inp / source["file"]).filter(
        (pl.col("endpoint") == source["endpoint"]) & (pl.col("method") == source["source_method"])
    )
    rows = frame.select(
        "compound_id",
        "endpoint",
        pl.lit(source["member"]).alias("member"),
        "prediction",
    )
    if rows.height != 750:
        raise ValueError(f"Missing supplemental prediction: {source['member']}")
    rename_rows.append(rows)

atomic = pl.concat([atomic, *rename_rows])
if atomic.select("compound_id", "endpoint", "member").n_unique() != atomic.height:
    raise ValueError("Duplicate atomic prediction keys")


def weights_for(endpoint: str, recipe: dict) -> dict[str, float]:
    members = recipe["components"]
    if recipe["method"] == "equal_all":
        return {member: 1.0 / len(members) for member in members}
    fits = pl.read_csv(inp / recipe["fits_file"]).filter(pl.col("method") == recipe["method"])
    if fits.height != 5 or sorted(fits["outer_fold"].to_list()) != list(range(5)):
        raise ValueError(f"Invalid weight fits for {endpoint}")
    values = np.asarray([json.loads(value) for value in fits["weights"]], dtype=float)
    if values.shape != (5, len(members)):
        raise ValueError(f"Weight dimension mismatch for {endpoint}")
    weights = values.mean(axis=0)
    weights = weights / weights.sum()
    return {
        member: float(weight)
        for member, weight in zip(members, weights, strict=True)
        if weight > 1e-12
    }


endpoint_predictions = {}
weight_summary = {}
atomic_rows = []
for endpoint, recipe in cfg["recipes"].items():
    endpoint_weights = weights_for(endpoint, recipe)
    weight_summary[endpoint] = endpoint_weights
    available = atomic.filter(pl.col("endpoint") == endpoint)
    missing = set(endpoint_weights) - set(available["member"])
    if missing:
        raise ValueError((endpoint, sorted(missing)))
    prediction = np.zeros(test.height, dtype=float)
    for member, weight in endpoint_weights.items():
        frame = available.filter(pl.col("member") == member).select("compound_id", "prediction")
        aligned = test.select("compound_id").join(
            frame, on="compound_id", validate="1:1", maintain_order="left"
        )
        if aligned["prediction"].null_count():
            raise ValueError(f"Incomplete prediction: {endpoint}/{member}")
        values = aligned["prediction"].to_numpy()
        prediction += weight * values
        atomic_rows.append(
            test.select("compound_id").with_columns(
                pl.lit(endpoint).alias("endpoint"),
                pl.lit(member).alias("member"),
                pl.lit(weight).alias("weight"),
                pl.Series("prediction", values),
            )
        )
    endpoint_predictions[endpoint] = prediction

candidate = test.select("SMILES", "Molecule_Name")
for endpoint in cfg["recipes"]:
    candidate = candidate.with_columns(pl.Series(endpoint, endpoint_predictions[endpoint]))
candidate = candidate.sort("Molecule_Name")
candidate.write_csv(out / "candidate-raw.csv")
candidate.write_parquet(out / "candidate-raw.parquet")
pl.concat(atomic_rows).write_parquet(out / "atomic-test-predictions.parquet")

base = pl.read_csv(inp / cfg["base_candidate"]).sort("Molecule_Name")
if not candidate["Molecule_Name"].equals(base["Molecule_Name"]):
    raise ValueError("Candidate/base row mismatch")
placed = candidate
placement = {}
comparison = []
for endpoint in cfg["recipes"]:
    raw = candidate[endpoint].to_numpy()
    old = base[endpoint].to_numpy()
    adjusted = (raw - raw.mean()) / raw.std() * old.std() + old.mean()
    placed = placed.with_columns(pl.Series(endpoint, adjusted))
    placement[endpoint] = {
        "source_mean": float(raw.mean()),
        "source_std": float(raw.std()),
        "target_mean": float(old.mean()),
        "target_std": float(old.std()),
    }
    comparison.append(
        {
            "endpoint": endpoint,
            "members": len(weight_summary[endpoint]),
            "raw_mean_absolute_change": float(np.abs(raw - old).mean()),
            "placed_mean_absolute_change": float(np.abs(adjusted - old).mean()),
            "raw_base_spearman": float(spearmanr(raw, old).statistic),
            "placed_base_spearman": float(spearmanr(adjusted, old).statistic),
        }
    )
placed.write_csv(out / "candidate-placed.csv")
placed.write_parquet(out / "candidate-placed.parquet")
pl.DataFrame(comparison).write_csv(out / "candidate-comparison.csv")

for filename in ["candidate-raw.csv", "candidate-placed.csv"]:
    report = validate(out / filename, "regression", inp / cfg["reference"])
    if not report["valid"]:
        raise ValueError(report)

(out / "summary.json").write_text(
    json.dumps(
        {
            "recipes": cfg["recipes"],
            "weights": weight_summary,
            "comparison_to_previous_submission": comparison,
            "placement": placement,
            "all_four_endpoints_changed": True,
            "test_labels_used": False,
            "ensemble_topology": "target-specific flat atomic members",
        },
        indent=2,
    )
    + "\n"
)
print(json.dumps(comparison), flush=True)
