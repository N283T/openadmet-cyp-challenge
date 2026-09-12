"""Audit the packaged target-wise candidate with compressed CYP3A4 Caruana."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.tracks import validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
keys = ["compound_id", "endpoint"]

packaged = pl.read_parquet(inp / cfg["packaged_long"]).sort(keys)
assert packaged.height == 3000 and not packaged.select(keys).is_duplicated().any()
assert validate(inp / cfg["packaged_wide"], "regression", inp / cfg["reference"])["valid"]

components = {
    item["name"]: pl.read_parquet(inp / item["test"]).sort(keys)
    for item in cfg["components"]
}
for frame in components.values():
    assert frame.height == 3000 and frame["prediction"].is_finite().all()

errors = []
for endpoint, selected in cfg["endpoint_components"].items():
    expected = np.mean(
        np.stack(
            [
                components[name].filter(pl.col("endpoint") == endpoint)["prediction"].to_numpy()
                for name in selected
            ]
        ),
        axis=0,
    )
    actual = packaged.filter(pl.col("endpoint") == endpoint)["prediction"].to_numpy()
    errors.append(float(np.max(np.abs(actual - expected))))
assert max(errors) < 1e-12

previous = pl.read_parquet(inp / cfg["previous_packaged_long"]).sort(keys)
for endpoint in cfg["unchanged_endpoints"]:
    old = previous.filter(pl.col("endpoint") == endpoint)["prediction"].to_numpy()
    new = packaged.filter(pl.col("endpoint") == endpoint)["prediction"].to_numpy()
    assert np.allclose(old, new, rtol=0, atol=0)
old_3a4 = previous.filter(pl.col("endpoint") == cfg["endpoint"])["prediction"].to_numpy()
new_3a4 = packaged.filter(pl.col("endpoint") == cfg["endpoint"])["prediction"].to_numpy()
assert not np.allclose(old_3a4, new_3a4, rtol=0, atol=1e-12)

reconstruction = pl.read_csv(inp / cfg["reconstruction"])
assert reconstruction["max_abs_error"].max() < 1e-10
pruning = pl.read_csv(inp / cfg["pruning_scores"])
selected = pruning.filter(
    (pl.col("strategy") == "crossfold_topk")
    & (pl.col("weighting") == "weighted")
    & (pl.col("k") == 2)
).row(0, named=True)
assert selected["delta_ci_high"] < 0

top_level = pl.read_csv(inp / cfg["top_level_scores"])
cv_scores = {}
for endpoint, count in cfg["top_level_counts"].items():
    if endpoint == cfg["endpoint"]:
        cv_scores[endpoint] = selected["five_component_ST_RAE"]
    else:
        cv_scores[endpoint] = top_level.filter(
            (pl.col("endpoint") == endpoint) & (pl.col("component_count") == count)
        )["ST_RAE"].item()

verification = {
    "status": "passed",
    "test_means_reconstructed": True,
    "max_test_reconstruction_error": max(errors),
    "unchanged_endpoints_exact": cfg["unchanged_endpoints"],
    "cyp3a4_original_five_component_ST_RAE": selected["five_component_ST_RAE"]
    - selected["delta_vs_original_five"],
    "cyp3a4_top2_five_component_ST_RAE": selected["five_component_ST_RAE"],
    "cyp3a4_delta": selected["delta_vs_original_five"],
    "cyp3a4_delta_ci": [selected["delta_ci_low"], selected["delta_ci_high"]],
    "crossfit_endpoint_ST_RAE": cv_scores,
    "crossfit_macro_ST_RAE": float(np.mean(list(cv_scores.values()))),
    "caruana_internal_models_before": 77,
    "caruana_internal_models_after": cfg["caruana_internal_models_after"],
}
(out / "verification.json").write_text(json.dumps(verification, indent=2))
print(json.dumps(verification, indent=2), flush=True)
