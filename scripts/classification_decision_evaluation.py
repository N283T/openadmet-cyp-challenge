"""Evaluate model and call-rate choices on held-out classification folds."""

import json
import os
from pathlib import Path

import polars as pl

from cyp_challenge.classification_evaluation import evaluate_rates, prepare_predictions

inputs = Path(os.environ["CYP_INPUT_DIR"])
outputs = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

predictions = pl.read_parquet(inputs / cfg["predictions"])
renames = {"method": "model", "prediction": "score"}
predictions = predictions.rename(
    {old: new for old, new in renames.items() if old in predictions.columns}
)
assignments = pl.read_parquet(inputs / cfg["assignments"])
data = prepare_predictions(predictions, assignments, cfg["models"], cfg["strategy"])
summary, curve, decisions, calls = evaluate_rates(
    data,
    cfg["models"],
    [float(rate) for rate in cfg["rates"]],
    float(cfg["target_prior"]),
    int(cfg["expected_pool_rows"]),
    int(cfg["bootstrap_replicates"]),
    int(cfg["bootstrap_seed"]),
)
summary["endpoint"] = cfg["endpoint"]
summary["scope"] = cfg["scope"]
summary["interpretation"] = (
    "This is a labelled-only retrospective rate audit. It cannot reproduce "
    "all-compound call selection or justify a submission."
    if not summary["deployment_decision_supported"]
    else (
        "All-compound outer-fold pool is present; model selection still needs "
        "prospective validation."
    )
)
curve.write_csv(outputs / "rate-by-fold.csv")
decisions.write_csv(outputs / "nested-rate-choices.csv")
calls.write_parquet(outputs / "nested-calls.parquet")
(outputs / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
