"""Audit D6 OOF decision curves and blind changes across call counts."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.classification_evaluation import (
    fixed_rate_calls,
    prepare_predictions,
    projected_mcc,
)

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
bank = pl.read_parquet(inp / cfg["oof"])
assignments = pl.read_parquet(inp / cfg["assignments"])
rows = []
for strategy in cfg["strategies"]:
    data = prepare_predictions(
        bank.filter(pl.col("split_id") == strategy),
        assignments,
        cfg["oof_models"],
        strategy,
    )
    parts = [
        data.filter(pl.col("model") == model).sort("compound_id") for model in cfg["oof_models"]
    ]
    folds = parts[0]["fold"].to_numpy()
    y = parts[0]["y_true"].to_numpy()
    scored = np.isfinite(y)
    if parts[0].height != 6146 or scored.sum() != 1497:
        raise ValueError("Incomplete full-pool OOF")
    for count in cfg["call_counts"]:
        rate = count / cfg["blind_pool_rows"]
        for model, part in zip(cfg["oof_models"], parts, strict=True):
            calls = fixed_rate_calls(part["score"].to_numpy(), folds, rate)
            for prior in cfg["target_priors"]:
                rows.append(
                    {
                        "strategy": strategy,
                        "model": model,
                        "target_prior": prior,
                        "blind_call_count": count,
                        "oof_call_count": int(calls.sum()),
                        "scored_calls": int(calls[scored].sum()),
                        "scored_true_positives": int((calls[scored] & (y[scored] == 1)).sum()),
                        "projected_mcc": projected_mcc(y[scored].astype(int), calls[scored], prior),
                    }
                )
curve = pl.DataFrame(rows)
curve.write_csv(out / "oof-call-count-curve.csv")

blind = pl.read_parquet(inp / cfg["blind_scores"])
parts = [
    blind.filter(pl.col("model") == model).sort("compound_id") for model in cfg["blind_models"]
]
if any(part.height != cfg["blind_pool_rows"] for part in parts) or not parts[0][
    "compound_id"
].equals(parts[1]["compound_id"]):
    raise ValueError("Blind score alignment failed")
scores = {
    cfg["oof_models"][0]: parts[0]["score"].to_numpy(),
    cfg["oof_models"][1]: (1 - cfg["ap_weight"]) * parts[0]["score"].to_numpy()
    + cfg["ap_weight"] * parts[1]["score"].to_numpy(),
}
base = pl.read_parquet(inp / cfg["base_submission"]).select(
    "Molecule_Name", pl.col(cfg["endpoint"]).alias("submitted_call")
)
aliases = (
    pl.read_parquet(inp / cfg["aliases"])
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
ordered = (
    parts[0]
    .select("compound_id")
    .join(aliases, on="compound_id", validate="1:1")
    .join(base, on="Molecule_Name", validate="1:1")
    .sort("compound_id")
)
if ordered.height != cfg["blind_pool_rows"]:
    raise ValueError("Blind IDs incomplete")
submitted = ordered["submitted_call"].to_numpy()
if submitted.sum() != cfg["historical_calls"]:
    raise ValueError("Historical call count changed")
blind_rows = []
for model, score in scores.items():
    ranking = np.argsort(-score, kind="stable")
    for count in cfg["call_counts"]:
        calls = np.zeros(len(score), dtype=bool)
        calls[ranking[:count]] = True
        blind_rows.append(
            {
                "model": model,
                "blind_call_count": count,
                "retained_from_sep15": int((calls & submitted).sum()),
                "removed_from_sep15": int((~calls & submitted).sum()),
                "added_to_sep15": int((calls & ~submitted).sum()),
            }
        )
    if model == cfg["oof_models"][0]:
        top_historical = np.zeros(len(score), dtype=bool)
        top_historical[ranking[: cfg["historical_calls"]]] = True
        if not np.array_equal(top_historical, submitted):
            raise ValueError("Blind control score does not reproduce submitted vector")
pl.DataFrame(blind_rows).write_csv(out / "blind-call-count-changes.csv")
