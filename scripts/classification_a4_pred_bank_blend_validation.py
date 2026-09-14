"""Cross-fit and grouped-bootstrap candidate A4 auxiliary rank blends."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    ranked = np.empty(len(values))
    for fold in np.unique(folds):
        mask = folds == fold
        ranked[mask] = rankdata(values[mask], method="average") / (mask.sum() + 1.0)
    return ranked


def fixed_calls(scores: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(positions))))
        result[positions[np.argsort(-scores[positions], kind="stable")[:count]]] = True
    return result


def expected_mcc(y: np.ndarray, selected: np.ndarray) -> float:
    positive = y == 1
    tpr = float(selected[positive].mean())
    fpr = float(selected[~positive].mean())
    prevalence = cfg["target_prior"]
    tp = prevalence * tpr
    fn = prevalence * (1 - tpr)
    fp = (1 - prevalence) * fpr
    tn = (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


base = read("caruana_oof").sort("compound_id")
new = read("new_oof").filter(pl.col("strategy") == "umap")
assignments = (
    read("assignments").filter(pl.col("strategy") == "umap").select("compound_id", "group")
)
base = base.join(assignments, on="compound_id", validate="1:1")
y = base["y_true"].to_numpy().astype(int)
folds = base["fold"].to_numpy().astype(int)
groups = base["group"].cast(pl.String).to_numpy()
base_rank = fold_rank(base[cfg["base_method"]].to_numpy(), folds)

candidates = {"baseline": base_rank}
for method in cfg["candidate_methods"]:
    aligned = (
        new.filter(pl.col("method") == method)
        .select("compound_id", "prediction")
        .join(base.select("compound_id"), on="compound_id", maintain_order="right")
    )
    candidate_rank = fold_rank(aligned["prediction"].to_numpy(), folds)
    for weight in cfg["weights"]:
        if weight > 0:
            candidates[f"{method}@{weight:g}"] = (1 - weight) * base_rank + weight * candidate_rank

crossfit = np.empty(len(y))
choice_rows = []
for heldout in np.unique(folds):
    selection = folds != heldout
    scored = []
    for name, prediction in candidates.items():
        selected = fixed_calls(prediction[selection], folds[selection])
        scored.append((expected_mcc(y[selection], selected), name))
    best_score, best_name = max(scored)
    mask = folds == heldout
    crossfit[mask] = candidates[best_name][mask]
    choice_rows.append(
        {"heldout_fold": int(heldout), "selected": best_name, "selection_MCC": best_score}
    )

baseline_selected = fixed_calls(base_rank, folds)
crossfit_selected = fixed_calls(crossfit, folds)
crossfit_summary = {
    "baseline_MCC_target_prior": expected_mcc(y, baseline_selected),
    "crossfit_MCC_target_prior": expected_mcc(y, crossfit_selected),
    "baseline_MCC_observed": float(matthews_corrcoef(y, baseline_selected)),
    "crossfit_MCC_observed": float(matthews_corrcoef(y, crossfit_selected)),
    "baseline_AUROC": float(roc_auc_score(y, base_rank)),
    "crossfit_AUROC": float(roc_auc_score(y, crossfit)),
    "baseline_AP": float(average_precision_score(y, base_rank)),
    "crossfit_AP": float(average_precision_score(y, crossfit)),
}

group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
rng = np.random.default_rng(cfg["bootstrap_seed"])
bootstrap_rows = []
for method in cfg["candidate_methods"]:
    for weight in cfg["bootstrap_weights"]:
        name = f"{method}@{weight:g}"
        prediction = candidates[name]
        selected = fixed_calls(prediction, folds)
        deltas = []
        for _ in range(cfg["bootstrap_replicates"]):
            sampled = rng.choice(group_values, size=len(group_values), replace=True)
            rows = np.concatenate([group_rows[group] for group in sampled])
            deltas.append(
                expected_mcc(y[rows], selected[rows])
                - expected_mcc(y[rows], baseline_selected[rows])
            )
        bootstrap_rows.append(
            {
                "candidate": name,
                "delta_MCC_target_prior": expected_mcc(y, selected)
                - expected_mcc(y, baseline_selected),
                "bootstrap_p025": float(np.quantile(deltas, 0.025)),
                "bootstrap_p50": float(np.quantile(deltas, 0.5)),
                "bootstrap_p975": float(np.quantile(deltas, 0.975)),
                "probability_delta_positive": float(np.mean(np.asarray(deltas) > 0)),
            }
        )

pl.DataFrame(choice_rows).write_csv(out / "crossfit-choices.csv")
bootstrap = pl.DataFrame(bootstrap_rows).sort("delta_MCC_target_prior", descending=True)
bootstrap.write_csv(out / "group-bootstrap.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "crossfit": crossfit_summary,
            "choices": choice_rows,
            "bootstrap_groups": len(group_values),
            "selection_note": (
                "Cross-fit selection excludes each evaluated fold. Bootstrap intervals remain "
                "conditional on the screened candidate grid."
            ),
        },
        indent=2,
    )
    + "\n"
)
print(json.dumps(crossfit_summary, indent=2), flush=True)
print(bootstrap, flush=True)
