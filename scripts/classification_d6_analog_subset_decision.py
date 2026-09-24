"""Select D6 calls within a pseudo analog-expansion test pool."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.classification_evaluation import prepare_predictions, projected_mcc

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
analog_ids = set(pl.read_parquet(inp / cfg["analogs"])["analog_id"].to_list())
assignments = pl.read_parquet(inp / cfg["assignments"])
rows = []
for strategy, files in cfg["predictions"].items():
    bank = pl.concat([pl.read_parquet(inp / name) for name in files]).filter(
        pl.col("split_id") == strategy
    )
    models = cfg["models"][strategy]
    aligned = prepare_predictions(bank, assignments, models, strategy)
    parts = [aligned.filter(pl.col("model") == name).sort("compound_id") for name in models]
    ids = parts[0]["compound_id"].to_list()
    if ids != parts[1]["compound_id"].to_list():
        raise ValueError("Model IDs differ")
    subset = np.asarray([compound_id in analog_ids for compound_id in ids])
    if int(subset.sum()) != len(analog_ids):
        raise ValueError("Pseudo analog IDs are incomplete")
    y = parts[0]["y_true"].to_numpy()
    scored = subset & np.isfinite(y)
    truth = y[scored].astype(int)
    count = int(round(cfg["rate"] * len(analog_ids)))
    calls = []
    for part in parts:
        score = part["score"].to_numpy()
        selected_rows = np.flatnonzero(subset)
        best = selected_rows[np.argsort(-score[selected_rows], kind="stable")[:count]]
        chosen = np.zeros(len(ids), dtype=bool)
        chosen[best] = True
        calls.append(chosen)
    control = projected_mcc(truth, calls[0][scored], cfg["target_prior"])
    ap = projected_mcc(truth, calls[1][scored], cfg["target_prior"])
    groups = parts[0]["group"].to_numpy()
    scored_rows = np.flatnonzero(scored)
    unique_groups = np.unique(groups[scored_rows])
    group_rows = {g: scored_rows[groups[scored_rows] == g] for g in unique_groups}
    rng = np.random.default_rng(cfg["bootstrap_seed"])
    deltas = []
    for _ in range(cfg["bootstrap_replicates"]):
        drawn = rng.choice(unique_groups, len(unique_groups), replace=True)
        sample_rows = np.concatenate([group_rows[g] for g in drawn])
        sample_truth = y[sample_rows].astype(int)
        delta = projected_mcc(
            sample_truth, calls[1][sample_rows], cfg["target_prior"]
        ) - projected_mcc(sample_truth, calls[0][sample_rows], cfg["target_prior"])
        if np.isfinite(delta):
            deltas.append(delta)
    rows.append(
        {
            "strategy": strategy,
            "analog_pool_rows": len(analog_ids),
            "calls_per_model": count,
            "scored_analogs": int(scored.sum()),
            "positive_analogs": int(truth.sum()),
            "control_calls_scored": int(calls[0][scored].sum()),
            "ap_calls_scored": int(calls[1][scored].sum()),
            "control_mcc": control,
            "ap_mcc": ap,
            "ap_minus_control": ap - control,
            "ci025": float(np.quantile(deltas, 0.025)),
            "ci975": float(np.quantile(deltas, 0.975)),
        }
    )
result = pl.DataFrame(rows)
result.write_csv(out / "analog-subset-decision.csv")
print(result, flush=True)
