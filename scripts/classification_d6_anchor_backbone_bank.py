"""Select a new-family-only flat D6 bank around anchored Jeremy and CheMeleon."""

import json
import os
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(key: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[key])


def rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = rank(values[mask])
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


sources = {
    "lgbm": (read("lgbm_oof"), read("lgbm_test")),
    "readout": (read("readout_oof"), read("readout_test")),
    "tabpfn": (read("tabpfn_oof"), read("tabpfn_test")),
}
reference = (
    sources["lgbm"][0]
    .filter(
        (pl.col("endpoint") == cfg["endpoint"])
        & (pl.col("split_id") == cfg["split_id"])
        & (pl.col("method") == "current_fixed_ratio")
    )
    .sort("compound_id")
)
ids = reference.select("compound_id")
folds = reference["fold"].to_numpy().astype(int)
y = reference["y_true"].to_numpy().astype(int)
groups = reference["group"].cast(pl.String).to_numpy()
control = fold_rank(reference["prediction"].to_numpy(), folds)

first_test = sources["lgbm"][1].filter(
    (pl.col("endpoint") == cfg["endpoint"])
    & (pl.col("split_id") == cfg["split_id"])
)
test_ids = first_test.select("compound_id").unique().sort("compound_id")


def align(source: str, method: str, test: bool = False) -> np.ndarray:
    table = sources[source][int(test)].filter(
        (pl.col("endpoint") == cfg["endpoint"])
        & (pl.col("split_id") == cfg["split_id"])
        & (pl.col("method") == method)
    )
    target = test_ids if test else ids
    joined = target.join(
        table.select("compound_id", "prediction"),
        on="compound_id",
        validate="1:1",
        maintain_order="left",
    )
    if joined["prediction"].null_count() or not joined["prediction"].is_finite().all():
        raise ValueError(f"Alignment failed for {source}/{method}, test={test}")
    return joined["prediction"].to_numpy()


names = [member["name"] for member in cfg["members"]]
oof_x = np.column_stack(
    [fold_rank(align(member["source"], member["method"]), folds) for member in cfg["members"]]
)
test_x = np.column_stack(
    [rank(align(member["source"], member["method"], test=True)) for member in cfg["members"]]
)


def subset_score(x: np.ndarray, labels: np.ndarray, split: np.ndarray, members: list[int]) -> float:
    return expected_mcc(labels, fixed_calls(x[:, members].mean(axis=1), split))


def forward_path(x: np.ndarray, labels: np.ndarray, split: np.ndarray) -> list[int]:
    selected: list[int] = []
    remaining = list(range(x.shape[1]))
    while remaining:
        trials = [
            subset_score(x, labels, split, selected + [candidate])
            for candidate in remaining
        ]
        best = max(range(len(remaining)), key=lambda i: (trials[i], -remaining[i]))
        selected.append(remaining.pop(best))
    return selected


allowed_sizes = [size for size in cfg["subset_sizes"] if size <= len(names)]
nested = np.empty(len(y), dtype=float)
selection_rows = []
outer_members = []
for held in np.unique(folds):
    outer_train = folds != held
    outer_valid = ~outer_train
    inner_folds = folds[outer_train]
    inner_predictions = {size: np.empty(outer_train.sum()) for size in allowed_sizes}
    for inner_held in np.unique(inner_folds):
        inner_train = inner_folds != inner_held
        inner_valid = ~inner_train
        path = forward_path(
            oof_x[outer_train][inner_train],
            y[outer_train][inner_train],
            inner_folds[inner_train],
        )
        for size in allowed_sizes:
            inner_predictions[size][inner_valid] = oof_x[outer_train][inner_valid][
                :, path[:size]
            ].mean(axis=1)
    inner_scores = {
        size: expected_mcc(y[outer_train], fixed_calls(values, inner_folds))
        for size, values in inner_predictions.items()
    }
    selected_size = max(allowed_sizes, key=lambda size: (inner_scores[size], -size))
    path = forward_path(oof_x[outer_train], y[outer_train], folds[outer_train])
    members = path[:selected_size]
    outer_members.append(members)
    nested[outer_valid] = oof_x[outer_valid][:, members].mean(axis=1)
    selection_rows.append(
        {
            "held_fold": int(held),
            "selected_size": selected_size,
            "selected_members": json.dumps([names[index] for index in members]),
            "inner_scores": json.dumps(inner_scores, sort_keys=True),
        }
    )

counts = Counter(index for members in outer_members for index in members)
consensus = [
    index for index, count in counts.items() if count >= cfg["consensus_min_folds"]
]
if not consensus:
    consensus = [counts.most_common(1)[0][0]]
consensus.sort(key=lambda index: (-counts[index], index))

full_path = forward_path(oof_x, y, folds)
full_scores = {
    size: subset_score(oof_x, y, folds, full_path[:size]) for size in allowed_sizes
}
full_size = max(allowed_sizes, key=lambda size: (full_scores[size], -size))
full_members = full_path[:full_size]
primary_index = names.index(cfg["primary_member"])

predictions = {
    "current_control": control,
    "new_primary": oof_x[:, primary_index],
    "new_equal_all": oof_x.mean(axis=1),
    "new_nested_equal_subset": nested,
    "new_consensus_subset_descriptive": oof_x[:, consensus].mean(axis=1),
    "new_full_forward_descriptive": oof_x[:, full_members].mean(axis=1),
}
test_predictions = {
    "new_primary": test_x[:, primary_index],
    "new_equal_all": test_x.mean(axis=1),
    "new_consensus_subset": test_x[:, consensus].mean(axis=1),
    "new_full_forward": test_x[:, full_members].mean(axis=1),
}

metric_rows, call_map, prediction_frames = [], {}, []
for method, values in predictions.items():
    calls = fixed_calls(values, folds)
    call_map[method] = calls
    metric_rows.append(
        {
            "method": method,
            "AUROC": float(roc_auc_score(y, values)),
            "AP": float(average_precision_score(y, values)),
            "MCC_observed": float(matthews_corrcoef(y, calls)),
            "MCC_target_prior": expected_mcc(y, calls),
            "positive_calls": int(calls.sum()),
            "swaps_from_control": int(np.count_nonzero(calls != fixed_calls(control, folds))),
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

test_frames = []
for method, values in test_predictions.items():
    calls = np.zeros(len(values), dtype=bool)
    count = int(round(cfg["call_rate"] * len(values)))
    calls[np.argsort(-values, kind="stable")[:count]] = True
    test_frames.append(
        test_ids.with_columns(
            pl.lit(cfg["endpoint"]).alias("endpoint"),
            pl.lit(cfg["split_id"]).alias("split_id"),
            pl.lit(-1).alias("fold"),
            pl.lit(method).alias("method"),
            pl.Series("prediction", values),
            pl.Series("call", calls),
        )
    )

rng = np.random.default_rng(cfg["bootstrap_seed"])
unique_groups = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in unique_groups}
bootstrap_rows = []
for candidate in ["new_primary", "new_equal_all", "new_nested_equal_subset"]:
    deltas = []
    for _ in range(cfg["bootstrap_replicates"]):
        sampled = rng.choice(unique_groups, len(unique_groups), replace=True)
        rows = np.concatenate([group_rows[group] for group in sampled])
        deltas.append(
            expected_mcc(y[rows], call_map[candidate][rows])
            - expected_mcc(y[rows], call_map["current_control"][rows])
        )
    values = np.asarray(deltas)
    bootstrap_rows.append(
        {
            "candidate": candidate,
            "control": "current_control",
            "mean_delta": float(values.mean()),
            "ci025": float(np.quantile(values, 0.025)),
            "ci975": float(np.quantile(values, 0.975)),
            "probability_positive": float((values > 0).mean()),
        }
    )

standalone_rows = []
for index, name in enumerate(names):
    calls = fixed_calls(oof_x[:, index], folds)
    standalone_rows.append(
        {
            "member": name,
            "AUROC": float(roc_auc_score(y, oof_x[:, index])),
            "AP": float(average_precision_score(y, oof_x[:, index])),
            "MCC_target_prior": expected_mcc(y, calls),
            "control_spearman": float(spearmanr(oof_x[:, index], control).statistic),
            "selection_frequency": counts[index],
        }
    )

pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True).write_csv(
    out / "metrics.csv"
)
pl.DataFrame(standalone_rows).sort("MCC_target_prior", descending=True).write_csv(
    out / "bank-members.csv"
)
pl.DataFrame(selection_rows).write_parquet(out / "selections.parquet")
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
pl.concat(prediction_frames).write_parquet(out / "oof.parquet")
pl.concat(test_frames).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "protocol": cfg["protocol"],
            "members": names,
            "consensus_members": [names[index] for index in consensus],
            "full_forward_members": [names[index] for index in full_members],
            "full_forward_scores": full_scores,
            "candidate_created": False,
            "uploaded": False,
        },
        indent=2,
    )
    + "\n"
)
print(pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True), flush=True)
