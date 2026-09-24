"""Evaluate a hypothetical assay-interval mask without changing call selection."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.classification_evaluation import projected_mcc

inputs = Path(os.environ["CYP_INPUT_DIR"])
outputs = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
decisions = pl.read_parquet(inputs / cfg["decisions"])
route = (
    pl.read_parquet(inputs / cfg["assay_route"])
    .filter(pl.col("endpoint") == cfg["endpoint"])
    .select("compound_id", "shift_lower", "shift_upper")
)
if route.select("compound_id").is_duplicated().any():
    raise ValueError("Duplicate assay-route compound")
data = decisions.join(route, on="compound_id", how="left", validate="m:1")
assignments = pl.read_parquet(inputs / cfg["assignments"]).select(
    pl.col("strategy"), pl.col("compound_id"), pl.col("group")
)
data = data.join(assignments, on=["strategy", "compound_id"], how="left", validate="1:1")
if data["shift_lower"].null_count() == data.height:
    raise ValueError("No assay intervals joined")

rows = []
for strategy in sorted(data["strategy"].unique().to_list()):
    part = data.filter(pl.col("strategy") == strategy).sort("compound_id")
    y = part["y_true"].to_numpy()
    scored = np.isfinite(y)
    confident = (part["shift_lower"].to_numpy() > cfg["shift_threshold"]) | (
        part["shift_upper"].to_numpy() < cfg["shift_threshold"]
    )
    calls = {name: part[f"{name}_call"].to_numpy() for name in ["control", "ap"]}
    if part["group"].null_count():
        raise ValueError("Grouped uncertainty requires group IDs")
    groups = part["group"].to_numpy()
    for mask_name, mask in [
        ("all_supplied_labels", scored),
        ("shift_interval_excludes_threshold", scored & confident),
    ]:
        included = np.flatnonzero(mask)
        truth = y[included].astype(int)
        unique_groups = np.unique(groups[included])
        group_rows = {g: included[groups[included] == g] for g in unique_groups}
        rng = np.random.default_rng(cfg["bootstrap_seed"])
        samples = []
        for _ in range(cfg["bootstrap_replicates"]):
            drawn = rng.choice(unique_groups, len(unique_groups), replace=True)
            subset = np.concatenate([group_rows[g] for g in drawn])
            sample_truth = y[subset].astype(int)
            delta = projected_mcc(sample_truth, calls["ap"][subset], cfg["prior"]) - projected_mcc(
                sample_truth, calls["control"][subset], cfg["prior"]
            )
            if np.isfinite(delta):
                samples.append(delta)
        if not samples:
            raise ValueError("No valid bootstrap samples")
        control = projected_mcc(truth, calls["control"][included], cfg["prior"])
        ap = projected_mcc(truth, calls["ap"][included], cfg["prior"])
        rows.append(
            {
                "strategy": strategy,
                "mask": mask_name,
                "scored_rows": len(included),
                "positives": int(truth.sum()),
                "control_calls_scored": int(calls["control"][included].sum()),
                "ap_calls_scored": int(calls["ap"][included].sum()),
                "control_mcc": control,
                "ap_mcc": ap,
                "ap_minus_control": ap - control,
                "ci025": float(np.quantile(samples, 0.025)),
                "ci975": float(np.quantile(samples, 0.975)),
            }
        )
result = pl.DataFrame(rows)
result.write_csv(outputs / "scoring-mask-sensitivity.csv")
print(result, flush=True)
