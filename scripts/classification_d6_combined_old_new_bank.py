"""Compare and select a combined legacy plus anchored-family D6 bank."""

import json
import os
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata, spearmanr
from sklearn.linear_model import LogisticRegression
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


def metrics(values: np.ndarray, labels: np.ndarray, splits: np.ndarray) -> dict[str, float]:
    calls = fixed_calls(values, splits)
    return {
        "AUROC": float(roc_auc_score(labels, values)),
        "AP": float(average_precision_score(labels, values)),
        "MCC_target_prior": expected_mcc(labels, calls),
    }


legacy = read("legacy_bank").filter(pl.col("endpoint") == cfg["endpoint"])
reference = legacy.filter(pl.col("model") == cfg["legacy_control"]).sort("compound_id")
ids = reference.select("compound_id")
folds = reference["fold"].to_numpy().astype(int)
y = reference["y_true"].to_numpy().astype(int)
control = fold_rank(reference["prediction"].to_numpy(), folds)

legacy_names = sorted(
    name for name in legacy["model"].unique().to_list() if name != cfg["legacy_control"]
)
names = [f"legacy::{name}" for name in legacy_names]
origins = ["legacy"] * len(legacy_names)
values = []
for name in legacy_names:
    table = legacy.filter(pl.col("model") == name).select("compound_id", "prediction")
    joined = ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    if joined["prediction"].null_count():
        raise ValueError(f"Legacy alignment failed for {name}")
    values.append(fold_rank(joined["prediction"].to_numpy(), folds))

new_sources = {
    "new_lgbm": read("new_lgbm_oof"),
    "new_readout": read("new_readout_oof"),
    "new_tabpfn": read("new_tabpfn_oof"),
}
for spec in cfg["new_members"]:
    table = new_sources[spec["source"]].filter(
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
        raise ValueError(f"New alignment failed for {spec['name']}")
    names.append(spec["name"])
    origins.append("new")
    values.append(fold_rank(joined["prediction"].to_numpy(), folds))

x = np.column_stack(values)
legacy_indices = [index for index, origin in enumerate(origins) if origin == "legacy"]
new_indices = [index for index, origin in enumerate(origins) if origin == "new"]


def score_subset(
    matrix: np.ndarray,
    labels: np.ndarray,
    splits: np.ndarray,
    members: list[int],
) -> float:
    return expected_mcc(labels, fixed_calls(matrix[:, members].mean(axis=1), splits))


def forward_path(matrix: np.ndarray, labels: np.ndarray, splits: np.ndarray) -> list[int]:
    selected: list[int] = []
    remaining = list(range(matrix.shape[1]))
    while remaining:
        scores = [
            score_subset(matrix, labels, splits, selected + [candidate])
            for candidate in remaining
        ]
        best = max(range(len(remaining)), key=lambda i: (scores[i], -remaining[i]))
        selected.append(remaining.pop(best))
    return selected


def quality_keep(member_metrics: dict[str, float]) -> bool:
    floors = cfg["quality_floor"]
    return any(member_metrics[key] >= floors[key] for key in floors)


def prune(
    matrix: np.ndarray,
    labels: np.ndarray,
    splits: np.ndarray,
    use_quality: bool,
) -> tuple[list[int], dict[int, dict[str, float]], dict[int, str]]:
    member_metrics = {
        index: metrics(matrix[:, index], labels, splits)
        for index in range(matrix.shape[1])
    }
    eligible = [
        index
        for index in range(matrix.shape[1])
        if not use_quality or quality_keep(member_metrics[index])
    ]
    order = sorted(
        eligible,
        key=lambda index: (
            -member_metrics[index]["MCC_target_prior"],
            -member_metrics[index]["AUROC"],
            -member_metrics[index]["AP"],
            names[index],
        ),
    )
    retained: list[int] = []
    reasons: dict[int, str] = {}
    for index in range(matrix.shape[1]):
        if index not in eligible:
            reasons[index] = "quality_floor"
    for index in order:
        duplicate = next(
            (
                prior
                for prior in retained
                if spearmanr(matrix[:, index], matrix[:, prior]).statistic
                >= cfg["correlation_ceiling"]
            ),
            None,
        )
        if duplicate is None:
            retained.append(index)
        else:
            reasons[index] = f"correlated_with::{names[duplicate]}"
    return retained, member_metrics, reasons


predictions = {
    "current_control": control,
    "legacy_equal_all": x[:, legacy_indices].mean(axis=1),
    "new_equal_all": x[:, new_indices].mean(axis=1),
    "combined_equal_all": x.mean(axis=1),
}
selection_rows = []
frequency: Counter[int] = Counter()
for mode, use_quality in [("correlation_pruned", False), ("quality_correlation_pruned", True)]:
    prediction = np.empty(len(y), dtype=float)
    for held in np.unique(folds):
        train = folds != held
        valid = ~train
        retained, member_metrics, reasons = prune(x[train], y[train], folds[train], use_quality)
        prediction[valid] = x[valid][:, retained].mean(axis=1)
        if mode == "quality_correlation_pruned":
            frequency.update(retained)
        selection_rows.append(
            {
                "method": mode,
                "held_fold": int(held),
                "selected_size": len(retained),
                "selected_members": json.dumps([names[index] for index in retained]),
                "removed": json.dumps(
                    {names[index]: reason for index, reason in reasons.items()}, sort_keys=True
                ),
                "training_metrics": json.dumps(
                    {names[index]: value for index, value in member_metrics.items()},
                    sort_keys=True,
                ),
            }
        )
    predictions[f"combined_{mode}"] = prediction

allowed_sizes = [size for size in cfg["subset_sizes"] if size <= x.shape[1]]
nested = np.empty(len(y), dtype=float)
for held in np.unique(folds):
    outer_train = folds != held
    outer_valid = ~outer_train
    inner_folds = folds[outer_train]
    inner_predictions = {size: np.empty(outer_train.sum()) for size in allowed_sizes}
    for inner_held in np.unique(inner_folds):
        inner_train = inner_folds != inner_held
        inner_valid = ~inner_train
        path = forward_path(
            x[outer_train][inner_train],
            y[outer_train][inner_train],
            inner_folds[inner_train],
        )
        for size in allowed_sizes:
            inner_predictions[size][inner_valid] = x[outer_train][inner_valid][
                :, path[:size]
            ].mean(axis=1)
    inner_scores = {
        size: expected_mcc(y[outer_train], fixed_calls(values, inner_folds))
        for size, values in inner_predictions.items()
    }
    selected_size = max(allowed_sizes, key=lambda size: (inner_scores[size], -size))
    path = forward_path(x[outer_train], y[outer_train], folds[outer_train])
    retained = path[:selected_size]
    nested[outer_valid] = x[outer_valid][:, retained].mean(axis=1)
    selection_rows.append(
        {
            "method": "nested_equal_subset",
            "held_fold": int(held),
            "selected_size": selected_size,
            "selected_members": json.dumps([names[index] for index in retained]),
            "removed": json.dumps({}),
            "training_metrics": json.dumps(inner_scores, sort_keys=True),
        }
    )
predictions["combined_nested_equal_subset"] = nested

meta_rows = []
for penalty in cfg["meta_penalties"]:
    prediction = np.empty(len(y), dtype=float)
    for held in np.unique(folds):
        outer_train = folds != held
        outer_valid = ~outer_train
        inner_folds = folds[outer_train]
        inner_predictions = {
            value: np.empty(outer_train.sum(), dtype=float) for value in cfg["meta_Cs"]
        }
        for inner_held in np.unique(inner_folds):
            inner_train = inner_folds != inner_held
            inner_valid = ~inner_train
            for value in cfg["meta_Cs"]:
                model = LogisticRegression(
                    C=value,
                    penalty=penalty,
                    solver="liblinear",
                    class_weight="balanced",
                    max_iter=3000,
                    random_state=cfg["bootstrap_seed"] + int(held),
                )
                model.fit(x[outer_train][inner_train], y[outer_train][inner_train])
                inner_predictions[value][inner_valid] = model.predict_proba(
                    x[outer_train][inner_valid]
                )[:, 1]
        inner_scores = {
            value: expected_mcc(y[outer_train], fixed_calls(values, inner_folds))
            for value, values in inner_predictions.items()
        }
        selected_c = max(cfg["meta_Cs"], key=lambda value: (inner_scores[value], -value))
        model = LogisticRegression(
            C=selected_c,
            penalty=penalty,
            solver="liblinear",
            class_weight="balanced",
            max_iter=3000,
            random_state=cfg["bootstrap_seed"] + int(held),
        )
        model.fit(x[outer_train], y[outer_train])
        prediction[outer_valid] = model.predict_proba(x[outer_valid])[:, 1]
        for index, coefficient in enumerate(model.coef_[0]):
            meta_rows.append(
                {
                    "penalty": penalty,
                    "held_fold": int(held),
                    "selected_C": selected_c,
                    "member": names[index],
                    "origin": origins[index],
                    "coefficient": float(coefficient),
                    "inner_scores": json.dumps(inner_scores, sort_keys=True),
                }
            )
    predictions[f"combined_logistic_{penalty}"] = prediction

metric_rows, calls_by_method, oof_frames = [], {}, []
control_calls = fixed_calls(control, folds)
for method, prediction in predictions.items():
    calls = fixed_calls(prediction, folds)
    calls_by_method[method] = calls
    scored = metrics(prediction, y, folds)
    metric_rows.append(
        {
            "method": method,
            **scored,
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

rng = np.random.default_rng(cfg["bootstrap_seed"])
group_source = new_sources["new_lgbm"].filter(
    (pl.col("endpoint") == cfg["endpoint"])
    & (pl.col("split_id") == cfg["split_id"])
    & (pl.col("method") == cfg["new_members"][0]["method"])
)
group_join = ids.join(
    group_source.select("compound_id", "group"),
    on="compound_id",
    validate="1:1",
    maintain_order="left",
)
groups = group_join["group"].cast(pl.String).to_numpy()
unique_groups = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in unique_groups}
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

full_retained, full_member_metrics, full_reasons = prune(x, y, folds, True)
member_rows = []
for index, name in enumerate(names):
    member_rows.append(
        {
            "member": name,
            "origin": origins[index],
            **full_member_metrics[index],
            "mean_correlation": float(
                np.delete(np.corrcoef(x, rowvar=False)[index], index).mean()
            ),
            "max_correlation": float(
                np.delete(np.corrcoef(x, rowvar=False)[index], index).max()
            ),
            "full_prune_status": full_reasons.get(index, "kept"),
            "outer_keep_frequency": frequency[index],
        }
    )

cross_rows = []
correlation = np.corrcoef(x, rowvar=False)
for old in legacy_indices:
    for new in new_indices:
        cross_rows.append(
            {
                "legacy_member": names[old],
                "new_member": names[new],
                "spearman": float(correlation[old, new]),
            }
        )

pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True).write_csv(
    out / "metrics.csv"
)
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
pl.DataFrame(member_rows).sort("MCC_target_prior", descending=True).write_csv(
    out / "members.csv"
)
pl.DataFrame(cross_rows).sort("spearman", descending=True).write_csv(
    out / "cross-family-correlations.csv"
)
pl.DataFrame(selection_rows).write_parquet(out / "selections.parquet")
pl.DataFrame(meta_rows).write_parquet(out / "meta-fits.parquet")
pl.concat(oof_frames).write_parquet(out / "oof.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "protocol": cfg["protocol"],
            "legacy_members": len(legacy_indices),
            "new_members": len(new_indices),
            "combined_members": len(names),
            "full_quality_correlation_members": [names[index] for index in full_retained],
            "full_removed": {names[index]: reason for index, reason in full_reasons.items()},
            "candidate_created": False,
            "uploaded": False,
        },
        indent=2,
    )
    + "\n"
)
print(pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True), flush=True)
