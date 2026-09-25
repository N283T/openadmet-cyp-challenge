"""Conditional parent-group bootstrap of fixed labelled TDI call changes."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
frame = pl.read_parquet(inp / cfg["scenario_oof"])
groups = pl.read_parquet(inp / cfg["assignments"]).select("compound_id", "group")
rng = np.random.default_rng(cfg["seed"])


def expected_mcc(y, calls, prior):
    positive = y == 1
    if not positive.any() or positive.all():
        return 0.0
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    tp, fn = prior * tpr, prior * (1 - tpr)
    fp, tn = (1 - prior) * fpr, (1 - prior) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


rows = []
for target, candidates in cfg["candidates"].items():
    table = frame.filter((pl.col("target") == target) & (pl.col("scenario") == "current"))
    table = table.join(groups, on="compound_id", validate="m:1").sort("compound_id")
    for candidate in candidates:
        part = frame.filter(
            (pl.col("target") == target) & (pl.col("scenario") == candidate)
        ).select("compound_id", pl.col("call").alias("candidate_call"))
        paired = table.join(part, on="compound_id", validate="1:1").sort("compound_id")
        assert paired.height == table.height
        y = paired["y_true"].to_numpy()
        control = paired["call"].to_numpy()
        alternative = paired["candidate_call"].to_numpy()
        group = paired["group"].to_numpy()
        unique = np.unique(group)
        indices = [np.flatnonzero(group == value) for value in unique]
        prior = cfg["target_prior"][target]
        whole = np.arange(len(y))

        def delta(index, y=y, control=control, alternative=alternative, prior=prior):
            return expected_mcc(y[index], alternative[index], prior) - expected_mcc(
                y[index], control[index], prior
            )

        draws = np.empty(cfg["replicates"], dtype=float)
        for index in range(len(draws)):
            sampled = rng.integers(len(unique), size=len(unique))
            selected = np.concatenate([indices[item] for item in sampled])
            draws[index] = delta(selected)
        rows.append({
            "target": target,
            "candidate": candidate,
            "delta_MCC": delta(whole),
            "q025": float(np.quantile(draws, 0.025)),
            "q975": float(np.quantile(draws, 0.975)),
            "fraction_positive": float(np.mean(draws > 0)),
            "call_swaps": int(np.sum(alternative & ~control)),
            "groups": len(unique),
        })
pl.DataFrame(rows).write_csv(out / "bootstrap.csv")
(out / "summary.json").write_text(json.dumps({"config": cfg}, indent=2) + "\n")
