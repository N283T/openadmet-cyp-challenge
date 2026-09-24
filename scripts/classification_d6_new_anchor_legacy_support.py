"""Add legacy supports to fixed new-family anchors with held-fold exclusion."""

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
definition = json.loads((inp / cfg["combined_definition"]).read_text())


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


legacy = read("legacy_bank").filter(pl.col("endpoint") == cfg["endpoint"])
reference = legacy.filter(pl.col("model") == cfg["legacy_control"]).sort("compound_id")
ids = reference.select("compound_id")
folds = reference["fold"].to_numpy().astype(int)
y = reference["y_true"].to_numpy().astype(int)
control = fold_rank(reference["prediction"].to_numpy(), folds)

legacy_names = sorted(
    name for name in legacy["model"].unique().to_list() if name != cfg["legacy_control"]
)
legacy_x = []
for name in legacy_names:
    table = legacy.filter(pl.col("model") == name).select("compound_id", "prediction")
    joined = ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    legacy_x.append(fold_rank(joined["prediction"].to_numpy(), folds))
legacy_x = np.column_stack(legacy_x)

sources = {
    "new_lgbm": read("new_lgbm_oof"),
    "new_readout": read("new_readout_oof"),
    "new_tabpfn": read("new_tabpfn_oof"),
}
new_names = [spec["name"] for spec in definition["new_members"]]
new_x = []
for spec in definition["new_members"]:
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
    new_x.append(fold_rank(joined["prediction"].to_numpy(), folds))
new_x = np.column_stack(new_x)
new_index = {name: index for index, name in enumerate(new_names)}


def combine(
    anchor_matrix: np.ndarray,
    support_matrix: np.ndarray,
    support: list[int],
    mode: str,
) -> np.ndarray:
    if not support:
        return anchor_matrix.mean(axis=1)
    if mode == "atomic":
        return np.column_stack([anchor_matrix, support_matrix[:, support]]).mean(axis=1)
    if mode == "family_block":
        return np.column_stack(
            [anchor_matrix.mean(axis=1), support_matrix[:, support]]
        ).mean(axis=1)
    raise ValueError(mode)


def forward_path(
    anchor_matrix: np.ndarray,
    support_matrix: np.ndarray,
    labels: np.ndarray,
    split: np.ndarray,
    mode: str,
) -> list[int]:
    selected: list[int] = []
    remaining = list(range(support_matrix.shape[1]))
    while remaining:
        scores = [
            expected_mcc(
                labels,
                fixed_calls(
                    combine(anchor_matrix, support_matrix, selected + [candidate], mode),
                    split,
                ),
            )
            for candidate in remaining
        ]
        best = max(range(len(remaining)), key=lambda index: (scores[index], -remaining[index]))
        selected.append(remaining.pop(best))
    return selected


predictions = {"current_control": control}
selection_rows = []
support_frequency: dict[str, Counter[int]] = {}
for anchor_name, members in cfg["anchors"].items():
    anchor_indices = (
        list(range(len(new_names)))
        if members == "all"
        else [new_index[name] for name in members]
    )
    anchor_x = new_x[:, anchor_indices]
    predictions[f"{anchor_name}::no_support"] = anchor_x.mean(axis=1)
    for mode in cfg["weighting_modes"]:
        method = f"{anchor_name}::{mode}"
        prediction = np.empty(len(y), dtype=float)
        support_frequency[method] = Counter()
        for held in np.unique(folds):
            outer_train = folds != held
            outer_valid = ~outer_train
            inner_folds = folds[outer_train]
            inner_predictions = {
                size: np.empty(outer_train.sum()) for size in cfg["support_sizes"]
            }
            for inner_held in np.unique(inner_folds):
                inner_train = inner_folds != inner_held
                inner_valid = ~inner_train
                path = forward_path(
                    anchor_x[outer_train][inner_train],
                    legacy_x[outer_train][inner_train],
                    y[outer_train][inner_train],
                    inner_folds[inner_train],
                    mode,
                )
                for size in cfg["support_sizes"]:
                    inner_predictions[size][inner_valid] = combine(
                        anchor_x[outer_train][inner_valid],
                        legacy_x[outer_train][inner_valid],
                        path[:size],
                        mode,
                    )
            inner_scores = {
                size: expected_mcc(y[outer_train], fixed_calls(values, inner_folds))
                for size, values in inner_predictions.items()
            }
            selected_size = max(
                cfg["support_sizes"], key=lambda size: (inner_scores[size], -size)
            )
            path = forward_path(
                anchor_x[outer_train],
                legacy_x[outer_train],
                y[outer_train],
                folds[outer_train],
                mode,
            )
            support = path[:selected_size]
            support_frequency[method].update(support)
            prediction[outer_valid] = combine(
                anchor_x[outer_valid], legacy_x[outer_valid], support, mode
            )
            selection_rows.append(
                {
                    "method": method,
                    "held_fold": int(held),
                    "selected_size": selected_size,
                    "selected_support": json.dumps(
                        [legacy_names[index] for index in support]
                    ),
                    "inner_scores": json.dumps(inner_scores, sort_keys=True),
                }
            )
        predictions[method] = prediction

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

group_source = sources["new_lgbm"].filter(
    (pl.col("endpoint") == cfg["endpoint"])
    & (pl.col("split_id") == cfg["split_id"])
    & (pl.col("method") == definition["new_members"][0]["method"])
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
for candidate in predictions:
    if candidate == "current_control":
        continue
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

frequency_rows = []
for method, counts in support_frequency.items():
    for index, name in enumerate(legacy_names):
        frequency_rows.append(
            {"method": method, "legacy_support": name, "outer_frequency": counts[index]}
        )

metrics_table = pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True)
metrics_table.write_csv(out / "metrics.csv")
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
pl.DataFrame(selection_rows).write_parquet(out / "selections.parquet")
pl.DataFrame(frequency_rows).sort(
    "method", "outer_frequency", descending=[False, True]
).write_csv(out / "support-frequency.csv")
pl.concat(oof_frames).write_parquet(out / "oof.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "protocol": cfg["protocol"],
            "candidate_created": False,
            "uploaded": False,
        },
        indent=2,
    )
    + "\n"
)
print(metrics_table, flush=True)
