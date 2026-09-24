"""Run held-fold-excluded Caruana on one flat D6 old/backbone pool."""

import json
import os
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(key: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[key])


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = rankdata(values[mask], method="average") / (mask.sum() + 1.0)
    return result


def fixed_calls(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    calls = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(rows))))
        calls[rows[np.argsort(-values[rows], kind="stable")[:count]]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray) -> float:
    positive = y == 1
    if not positive.any() or positive.all():
        return 0.0
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    prevalence = cfg["target_prior"]
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def score(y: np.ndarray, values: np.ndarray, folds: np.ndarray, objective: str) -> float:
    if objective == "mcc":
        return expected_mcc(y, fixed_calls(values, folds))
    if objective == "auroc":
        return float(roc_auc_score(y, values))
    if objective == "ap":
        return float(average_precision_score(y, values))
    raise ValueError(objective)


legacy = read("legacy_bank").filter(pl.col("endpoint") == cfg["endpoint"])
reference = legacy.filter(pl.col("model") == cfg["control"]).sort("compound_id")
ids = reference.select("compound_id")
folds = reference["fold"].to_numpy().astype(int)
y = reference["y_true"].to_numpy().astype(int)
control = fold_rank(reference["prediction"].to_numpy(), folds)

names = []
values = []
for name in sorted(
    name for name in legacy["model"].unique().to_list() if name != cfg["control"]
):
    table = legacy.filter(pl.col("model") == name).select("compound_id", "prediction")
    joined = ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    names.append(name)
    values.append(fold_rank(joined["prediction"].to_numpy(), folds))

sources = {
    "lgbm": read("lgbm_oof"),
    "readout": read("readout_oof"),
    "tabpfn": read("tabpfn_oof"),
}
for spec in cfg["jeremy_chemeleon_members"]:
    table = sources[spec["source"]].filter(
        (pl.col("endpoint") == cfg["endpoint"])
        & (pl.col("split_id") == cfg["split_id"])
        & (pl.col("method") == spec["method"])
    )
    joined = ids.join(
        table.select("compound_id", "prediction"),
        on="compound_id",
        validate="1:1",
        maintain_order="left",
    )
    if joined["prediction"].null_count():
        raise ValueError(f"Alignment failed for {spec['name']}")
    names.append(spec["name"])
    values.append(fold_rank(joined["prediction"].to_numpy(), folds))
x = np.column_stack(values)


def fit_caruana(
    matrix: np.ndarray,
    labels: np.ndarray,
    split: np.ndarray,
    recipe: str,
) -> tuple[np.ndarray, list[dict[str, object]]]:
    if recipe == "equal_all":
        return np.full(matrix.shape[1], 1 / matrix.shape[1]), []
    parts = recipe.split("_")
    objective = parts[1]
    cap = cfg["sparse_cap"] if recipe.endswith("sparse6") else matrix.shape[1]
    counts = np.zeros(matrix.shape[1], dtype=int)
    total = np.zeros(len(labels), dtype=float)
    best_score = -np.inf
    best_weights = None
    path = []
    for step in range(1, cfg["steps"] + 1):
        candidates = [
            index
            for index in range(matrix.shape[1])
            if counts[index] or np.count_nonzero(counts) < cap
        ]
        candidate_scores = [
            score(labels, (total + matrix[:, index]) / step, split, objective)
            for index in candidates
        ]
        position = int(np.argmax(candidate_scores))
        choice = candidates[position]
        total += matrix[:, choice]
        counts[choice] += 1
        value = candidate_scores[position]
        path.append(
            {"step": step, "added": names[choice], "objective": float(value)}
        )
        if value > best_score + 1e-12:
            best_score = value
            best_weights = counts.astype(float) / step
    if best_weights is None:
        raise RuntimeError("Caruana failed to select a member")
    return best_weights, path


predictions = {"current_control": control}
selection_rows = []
selection_frequency: dict[str, Counter[int]] = {}
for recipe in cfg["recipes"]:
    prediction = np.empty(len(y), dtype=float)
    selection_frequency[recipe] = Counter()
    for held in np.unique(folds):
        train = folds != held
        valid = ~train
        weights, path = fit_caruana(x[train], y[train], folds[train], recipe)
        prediction[valid] = x[valid] @ weights
        selected = np.flatnonzero(weights)
        selection_frequency[recipe].update(selected.tolist())
        selection_rows.append(
            {
                "recipe": recipe,
                "held_fold": int(held),
                "selected_members": json.dumps([names[index] for index in selected]),
                "weights": json.dumps(
                    {names[index]: float(weights[index]) for index in selected},
                    sort_keys=True,
                ),
                "path": json.dumps(path),
            }
        )
    predictions[recipe] = prediction

full_rows = []
for recipe in cfg["recipes"]:
    weights, path = fit_caruana(x, y, folds, recipe)
    for index, weight in enumerate(weights):
        full_rows.append(
            {
                "recipe": recipe,
                "member": names[index],
                "weight": float(weight),
                "outer_selection_frequency": selection_frequency[recipe][index],
                "full_selected": bool(weight > 0),
                "full_path": json.dumps(path) if index == 0 else "",
            }
        )

metric_rows, calls_by_method, oof_frames = [], {}, []
control_calls = fixed_calls(control, folds)
for method, prediction in predictions.items():
    calls = fixed_calls(prediction, folds)
    calls_by_method[method] = calls
    metric_rows.append(
        {
            "method": method,
            "AUROC": float(roc_auc_score(y, prediction)),
            "AP": float(average_precision_score(y, prediction)),
            "MCC_target_prior": expected_mcc(y, calls),
            "MCC_observed": float(matthews_corrcoef(y, calls)),
            "positive_calls": int(calls.sum()),
            "swaps_from_control": int(np.count_nonzero(calls != control_calls)),
        }
    )
    oof_frames.append(
        ids.with_columns(
            pl.Series("fold", folds),
            pl.Series("y_true", y),
            pl.lit(method).alias("method"),
            pl.Series("prediction", prediction),
            pl.Series("call", calls),
        )
    )

group_source = sources["lgbm"].filter(
    (pl.col("endpoint") == cfg["endpoint"])
    & (pl.col("split_id") == cfg["split_id"])
    & (pl.col("method") == cfg["jeremy_chemeleon_members"][0]["method"])
)
groups = ids.join(
    group_source.select("compound_id", "group"),
    on="compound_id",
    validate="1:1",
    maintain_order="left",
)["group"].cast(pl.String).to_numpy()
unique_groups = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in unique_groups}
rng = np.random.default_rng(cfg["bootstrap_seed"])
bootstrap_rows = []
for candidate in cfg["recipes"]:
    deltas = []
    for _ in range(cfg["bootstrap_replicates"]):
        sampled = rng.choice(unique_groups, len(unique_groups), replace=True)
        rows = np.concatenate([group_rows[group] for group in sampled])
        deltas.append(
            expected_mcc(y[rows], calls_by_method[candidate][rows])
            - expected_mcc(y[rows], control_calls[rows])
        )
    delta = np.asarray(deltas)
    bootstrap_rows.append(
        {
            "candidate": candidate,
            "control": "current_control",
            "mean_delta": float(delta.mean()),
            "ci025": float(np.quantile(delta, 0.025)),
            "ci975": float(np.quantile(delta, 0.975)),
            "probability_positive": float((delta > 0).mean()),
        }
    )

metrics_table = pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True)
metrics_table.write_csv(out / "metrics.csv")
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
pl.DataFrame(selection_rows).write_parquet(out / "outer-selections.parquet")
pl.DataFrame(full_rows).write_parquet(out / "full-selection.parquet")
pl.concat(oof_frames).write_parquet(out / "oof.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "protocol": cfg["protocol"],
            "pool_size": len(names),
            "members": names,
            "candidate_created": False,
            "uploaded": False,
        },
        indent=2,
    )
    + "\n"
)
print(metrics_table, flush=True)
