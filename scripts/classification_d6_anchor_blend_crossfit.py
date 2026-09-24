"""Cross-fit the fixed anchor+pub030 blend-weight screen."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
table = pl.read_parquet(inp / cfg["blend_oof"])
models = [cfg["base_model"]] + sorted(
    [
        model
        for model in table["model"].unique().to_list()
        if model.startswith(cfg["candidate_prefix"])
    ],
    key=lambda value: float(value.rsplit("w", 1)[1]),
)
if len(models) != 6:
    raise ValueError(f"Expected control plus five weights, found {models}")

reference = table.filter(pl.col("model") == cfg["base_model"]).sort("compound_id")
ids = reference["compound_id"].to_list()
folds = reference["fold"].to_numpy().astype(int)
groups = reference["group"].to_numpy()
y = reference["y_true"].to_numpy().astype(int)
predictions = {}
for model in models:
    model_table = table.filter(pl.col("model") == model).sort("compound_id")
    if model_table["compound_id"].to_list() != ids:
        raise ValueError(f"Identity mismatch for {model}")
    predictions[model] = model_table["prediction"].to_numpy()


def fixed_calls(values: np.ndarray, selected_rows: np.ndarray | None = None) -> np.ndarray:
    calls = np.zeros(len(values), dtype=bool)
    active_folds = np.unique(folds) if selected_rows is None else np.unique(folds[selected_rows])
    for fold in active_folds:
        rows = np.flatnonzero(folds == fold)
        if selected_rows is not None:
            rows = np.intersect1d(rows, selected_rows, assume_unique=True)
        count = max(1, int(round(cfg["call_rate"] * len(rows))))
        calls[rows[np.argsort(-values[rows], kind="stable")[:count]]] = True
    return calls


def expected_mcc(labels: np.ndarray, calls: np.ndarray) -> float:
    positive = labels == 1
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    prevalence = cfg["target_prior"]
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


crossfit_prediction = np.empty(len(y), dtype=float)
crossfit_calls = np.zeros(len(y), dtype=bool)
selection_rows = []
for held_fold in np.unique(folds):
    train_rows = np.flatnonzero(folds != held_fold)
    scores = []
    for model in models:
        calls = fixed_calls(predictions[model], train_rows)
        score = expected_mcc(y[train_rows], calls[train_rows])
        scores.append((score, -models.index(model), model))
    selected_score, _, selected_model = max(scores)
    held_rows = np.flatnonzero(folds == held_fold)
    crossfit_prediction[held_rows] = predictions[selected_model][held_rows]
    selected_calls = fixed_calls(predictions[selected_model])
    crossfit_calls[held_rows] = selected_calls[held_rows]
    selection_rows.append(
        {
            "held_fold": int(held_fold),
            "selected_model": selected_model,
            "training_MCC_target_prior": float(selected_score),
        }
    )

full_scores = []
for model in models:
    calls = fixed_calls(predictions[model])
    full_scores.append((expected_mcc(y, calls), -models.index(model), model))
full_score, _, full_model = max(full_scores)
selection_rows.append(
    {
        "held_fold": -1,
        "selected_model": full_model,
        "training_MCC_target_prior": float(full_score),
    }
)

base_prediction = predictions[cfg["base_model"]]
base_calls = fixed_calls(base_prediction)
rows = []
for model, prediction, calls in [
    (cfg["base_model"], base_prediction, base_calls),
    ("crossfit_weight_selection", crossfit_prediction, crossfit_calls),
]:
    rows.append(
        {
            "model": model,
            "AUROC": float(roc_auc_score(y, prediction)),
            "AP": float(average_precision_score(y, prediction)),
            "MCC_observed": float(matthews_corrcoef(y, calls)),
            "MCC_target_prior": expected_mcc(y, calls),
            "positive_calls": int(calls.sum()),
        }
    )

group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
rng = np.random.default_rng(cfg["bootstrap_seed"])
deltas = []
for _ in range(cfg["bootstrap_replicates"]):
    sampled = rng.choice(group_values, size=len(group_values), replace=True)
    sampled_rows = np.concatenate([group_rows[group] for group in sampled])
    deltas.append(
        expected_mcc(y[sampled_rows], crossfit_calls[sampled_rows])
        - expected_mcc(y[sampled_rows], base_calls[sampled_rows])
    )
delta = np.asarray(deltas)
bootstrap = {
    "model": "crossfit_weight_selection",
    "reference": cfg["base_model"],
    "mean_MCC_delta": float(delta.mean()),
    "ci025": float(np.quantile(delta, 0.025)),
    "ci975": float(np.quantile(delta, 0.975)),
    "probability_positive": float((delta > 0).mean()),
}

pl.DataFrame(rows).write_csv(out / "metrics.csv")
pl.DataFrame(selection_rows).write_csv(out / "selections.csv")
pl.DataFrame([bootstrap]).write_csv(out / "bootstrap.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "selection": cfg["selection"],
            "full_data_choice": full_model,
            "candidate_created": False,
            "uploaded": False,
        },
        indent=2,
    )
    + "\n"
)
print(pl.DataFrame(rows), flush=True)
print(pl.DataFrame(selection_rows), flush=True)
