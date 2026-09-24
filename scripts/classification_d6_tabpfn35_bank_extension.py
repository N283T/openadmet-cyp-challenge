"""Add TabPFN 3.5 to the D6 bank without manual blend weights."""

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


def read(key: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[key])


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = unit_rank(values[mask])
    return result


def fixed_calls(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    calls = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(rows))))
        order = np.argsort(-values[rows], kind="stable")
        calls[rows[order[:count]]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray) -> float:
    positive = y == 1
    if not positive.any() or positive.all():
        return 0.0
    prevalence = cfg["target_prior"]
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def score_subset(x: np.ndarray, y: np.ndarray, folds: np.ndarray, members: list[int]) -> float:
    prediction = x[:, members].mean(axis=1)
    return expected_mcc(y, fixed_calls(prediction, folds))


def forward_path(x: np.ndarray, y: np.ndarray, folds: np.ndarray) -> list[int]:
    selected: list[int] = []
    remaining = list(range(x.shape[1]))
    while remaining:
        trials = [score_subset(x, y, folds, selected + [candidate]) for candidate in remaining]
        best_position = int(np.argmax(trials))
        selected.append(remaining.pop(best_position))
    return selected


def nested_equal_subset(
    x: np.ndarray, y: np.ndarray, folds: np.ndarray, names: list[str], bank_name: str
) -> tuple[np.ndarray, list[dict[str, object]]]:
    prediction = np.empty(len(y), dtype=float)
    records: list[dict[str, object]] = []
    allowed_sizes = [size for size in cfg["subset_sizes"] if size <= x.shape[1]]
    for held in np.unique(folds):
        outer_train = folds != held
        outer_valid = folds == held
        inner_folds = folds[outer_train]
        inner_prediction = {
            size: np.empty(outer_train.sum(), dtype=float) for size in allowed_sizes
        }
        for inner_held in np.unique(inner_folds):
            inner_train = inner_folds != inner_held
            inner_valid = inner_folds == inner_held
            path = forward_path(
                x[outer_train][inner_train],
                y[outer_train][inner_train],
                inner_folds[inner_train],
            )
            for size in allowed_sizes:
                members = path[:size]
                inner_prediction[size][inner_valid] = x[outer_train][inner_valid][
                    :, members
                ].mean(axis=1)
        inner_scores = {
            size: expected_mcc(y[outer_train], fixed_calls(values, inner_folds))
            for size, values in inner_prediction.items()
        }
        selected_size = max(allowed_sizes, key=lambda size: (inner_scores[size], -size))
        path = forward_path(x[outer_train], y[outer_train], folds[outer_train])
        members = path[:selected_size]
        prediction[outer_valid] = x[outer_valid][:, members].mean(axis=1)
        records.append(
            {
                "bank": bank_name,
                "held_fold": int(held),
                "selected_size": selected_size,
                "selected_members": json.dumps([names[index] for index in members]),
                "inner_scores": json.dumps(inner_scores, sort_keys=True),
            }
        )
    return prediction, records


endpoint = cfg["endpoint"]
legacy_long = read("legacy_bank").filter(pl.col("endpoint") == endpoint)
base = (
    read("deployed_oof")
    .filter(pl.col("model") == cfg["deployed_model"])
    .sort("compound_id")
)
ids = base.select("compound_id")
folds = base["fold"].to_numpy().astype(int)
y = base["y_true"].to_numpy().astype(int)
groups = base["group"].cast(pl.String).to_numpy()

legacy_names = sorted(legacy_long["model"].unique().to_list())
legacy_names.remove(cfg["replace_legacy_model"])
legacy_names = [cfg["deployed_model"], *legacy_names]
legacy_values = [fold_rank(base["prediction"].to_numpy(), folds)]
for name in legacy_names[1:]:
    table = legacy_long.filter(pl.col("model") == name).select(
        "compound_id", "fold", "y_true", "prediction"
    )
    joined = ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    if (
        joined.height != ids.height
        or joined["prediction"].null_count()
        or not np.array_equal(joined["fold"].to_numpy(), folds)
        or not np.array_equal(joined["y_true"].to_numpy(), y)
    ):
        raise ValueError(f"Legacy bank alignment failed for {name}")
    legacy_values.append(fold_rank(joined["prediction"].to_numpy(), folds))
legacy_x = np.column_stack(legacy_values)

new_names: list[str] = []
new_values: list[np.ndarray] = []
tabpfn = read("tabpfn35_oof").filter(
    (pl.col("endpoint") == endpoint) & (pl.col("split_id") == "umap")
)
for method in cfg["new_methods"]:
    name = f"tabpfn35fast::{method}"
    table = tabpfn.filter(pl.col("method") == method).select(
        "compound_id", "fold", "y_true", "prediction"
    )
    joined = ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    if (
        joined.height != ids.height
        or joined["prediction"].null_count()
        or not np.array_equal(joined["fold"].to_numpy(), folds)
        or not np.array_equal(joined["y_true"].to_numpy(), y)
    ):
        raise ValueError(f"TabPFN bank alignment failed for {method}")
    new_names.append(name)
    new_values.append(fold_rank(joined["prediction"].to_numpy(), folds))
extended_names = legacy_names + new_names
extended_x = np.column_stack([legacy_x, *new_values])

predictions = {
    "deployed_mean_block": legacy_x[:, 0],
    "legacy_equal_mean_all": legacy_x.mean(axis=1),
    "extended_equal_mean_all": extended_x.mean(axis=1),
}
selection_rows: list[dict[str, object]] = []
for bank_name, names, x in [
    ("legacy", legacy_names, legacy_x),
    ("extended", extended_names, extended_x),
]:
    values, records = nested_equal_subset(x, y, folds, names, bank_name)
    predictions[f"{bank_name}_nested_equal_subset"] = values
    selection_rows.extend(records)

metric_rows = []
call_by_method = {}
prediction_frames = []
for method, values in predictions.items():
    calls = fixed_calls(values, folds)
    call_by_method[method] = calls
    metric_rows.append(
        {
            "method": method,
            "AUROC": float(roc_auc_score(y, values)),
            "AP": float(average_precision_score(y, values)),
            "MCC_observed": float(matthews_corrcoef(y, calls)),
            "MCC_target_prior": expected_mcc(y, calls),
            "positive_calls": int(calls.sum()),
        }
    )
    prediction_frames.append(
        ids.with_columns(
            pl.Series("fold", folds),
            pl.Series("group", groups),
            pl.Series("y_true", y),
            pl.lit(method).alias("method"),
            pl.Series("prediction", values),
            pl.Series("call", calls),
        )
    )

rng = np.random.default_rng(cfg["bootstrap_seed"])
group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
contrasts = [
    ("extended_equal_mean_all", "legacy_equal_mean_all"),
    ("extended_nested_equal_subset", "legacy_nested_equal_subset"),
    ("extended_nested_equal_subset", "deployed_mean_block"),
]
bootstrap_rows = []
for candidate, control in contrasts:
    deltas = []
    for _ in range(cfg["bootstrap_replicates"]):
        sampled = rng.choice(group_values, size=len(group_values), replace=True)
        rows = np.concatenate([group_rows[group] for group in sampled])
        deltas.append(
            expected_mcc(y[rows], call_by_method[candidate][rows])
            - expected_mcc(y[rows], call_by_method[control][rows])
        )
    delta = np.asarray(deltas)
    bootstrap_rows.append(
        {
            "candidate": candidate,
            "control": control,
            "mean_MCC_delta": float(delta.mean()),
            "ci025": float(np.quantile(delta, 0.025)),
            "ci975": float(np.quantile(delta, 0.975)),
            "probability_positive": float((delta > 0).mean()),
        }
    )

metrics = pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(selection_rows).write_parquet(out / "selections.parquet")
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
pl.concat(prediction_frames).write_parquet(out / "oof.parquet")
pl.DataFrame(
    {
        "bank": ["legacy"] * len(legacy_names) + ["extended_new"] * len(new_names),
        "model": legacy_names + new_names,
    }
).write_csv(out / "bank-members.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "protocol": cfg["protocol"],
            "legacy_experts": len(legacy_names),
            "extended_experts": len(extended_names),
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
print(pl.DataFrame(bootstrap_rows), flush=True)
