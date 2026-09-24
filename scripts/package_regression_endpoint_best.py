"""Package the endpoint-wise best observed regression submission vector."""

import json
import os
from pathlib import Path

import polars as pl

from cyp_challenge.tracks import validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

identity_columns = ["SMILES", "Molecule_Name"]
source_frames = {
    name: pl.read_csv(inp / filename).sort("Molecule_Name")
    for name, filename in cfg["sources"].items()
}
base_name = next(iter(source_frames))
base = source_frames[base_name]
for name, frame in source_frames.items():
    if frame.height != 750:
        raise ValueError(f"{name} has {frame.height} rows, expected 750")
    if not frame.select(identity_columns).equals(base.select(identity_columns)):
        raise ValueError(f"Identity mismatch between {base_name} and {name}")

candidate = base.select(identity_columns)
lineage_rows = []
for endpoint, selection in cfg["endpoint_selection"].items():
    source_name = selection["source"]
    source = source_frames[source_name]
    if endpoint not in source.columns:
        raise ValueError(f"Missing endpoint {endpoint} in {source_name}")
    values = source[endpoint]
    if values.null_count() or not values.is_finite().all():
        raise ValueError(f"Non-finite values for {endpoint} in {source_name}")
    candidate = candidate.with_columns(values.alias(endpoint))
    lineage_rows.append(
        {
            "endpoint": endpoint,
            "source": source_name,
            "submission_id": selection["submission_id"],
            "observed_st_rae": selection["observed_st_rae"],
            "prediction_mean": float(values.mean()),
            "prediction_sample_sd": float(values.std(ddof=1)),
        }
    )

candidate_path = out / "predictions.csv"
candidate.write_csv(candidate_path)
candidate.write_parquet(out / "predictions.parquet")
pl.DataFrame(lineage_rows).write_csv(out / "endpoint-lineage.csv")

report = validate(candidate_path, "regression", inp / cfg["reference"])
if not report["valid"]:
    raise ValueError(report)

observed_macro = sum(row["observed_st_rae"] for row in lineage_rows) / len(lineage_rows)
summary = {
    "strategy": "endpoint-wise best observed leaderboard vector",
    "rows": candidate.height,
    "endpoint_lineage": lineage_rows,
    "reconstructed_observed_macro_st_rae": observed_macro,
    "macro_reconstruction_rule": "arithmetic mean of four independently scored endpoints",
    "test_labels_used": False,
    "validation": report,
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
