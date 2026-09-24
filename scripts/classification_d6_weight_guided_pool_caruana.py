"""Adjust one D6 pool from inner-fold unrestricted Caruana weights."""

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
pool_cfg = json.loads((inp / cfg["pool_definition"]).read_text())


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
    if objective == "ap":
        return float(average_precision_score(y, values))
    if objective == "auroc":
        return float(roc_auc_score(y, values))
    raise ValueError(objective)


legacy = read("legacy_bank").filter(pl.col("endpoint") == cfg["endpoint"])
reference = legacy.filter(pl.col("model") == cfg["control"]).sort("compound_id")
ids = reference.select("compound_id")
folds = reference["fold"].to_numpy().astype(int)
y = reference["y_true"].to_numpy().astype(int)
control = fold_rank(reference["prediction"].to_numpy(), folds)

names: list[str] = []
values: list[np.ndarray] = []
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
for spec in pool_cfg["jeremy_chemeleon_members"]:
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
    objective: str,
    candidates: np.ndarray,
) -> np.ndarray:
    counts = np.zeros(matrix.shape[1], dtype=int)
    total = np.zeros(len(labels), dtype=float)
    best_score = -np.inf
    best_weights = None
    for step in range(1, cfg["steps"] + 1):
        candidate_scores = np.asarray(
            [
                score(labels, (total + matrix[:, index]) / step, split, objective)
                for index in candidates
            ]
        )
        best = np.flatnonzero(candidate_scores >= candidate_scores.max() - 1e-12)
        if objective == "mcc" and len(best) > 1:
            secondary = np.asarray(
                [
                    average_precision_score(
                        labels, (total + matrix[:, candidates[pos]]) / step
                    )
                    for pos in best
                ]
            )
            position = int(best[int(np.argmax(secondary))])
        else:
            position = int(best[0])
        choice = int(candidates[position])
        total += matrix[:, choice]
        counts[choice] += 1
        value = float(candidate_scores[position])
        if value > best_score + 1e-12:
            best_score = value
            best_weights = counts.astype(float) / step
    if best_weights is None:
        raise RuntimeError("Caruana failed")
    return best_weights


def adjust_pool(
    policy: dict[str, object],
    inner_weights: np.ndarray,
    train_matrix: np.ndarray,
) -> np.ndarray:
    frequency = (inner_weights > 0).sum(axis=0)
    mean_weight = inner_weights.mean(axis=0)
    if policy["name"] == "all":
        selected = np.arange(train_matrix.shape[1])
    elif "minimum_inner_frequency" in policy:
        selected = np.flatnonzero(frequency >= int(policy["minimum_inner_frequency"]))
    else:
        order = np.argsort(-mean_weight, kind="stable")
        positive = order[mean_weight[order] > 0]
        target = float(policy["cumulative_mean_weight"])
        cumulative = np.cumsum(mean_weight[positive])
        count = int(np.searchsorted(cumulative, target, side="left") + 1)
        selected = positive[:count]
    if not len(selected):
        selected = np.asarray([int(np.argmax(mean_weight))])
    threshold = policy.get("max_abs_correlation")
    if threshold is not None and len(selected) > 1:
        ordered = selected[np.argsort(-mean_weight[selected], kind="stable")]
        kept: list[int] = []
        for index in ordered:
            if not kept:
                kept.append(int(index))
                continue
            correlations = np.corrcoef(
                train_matrix[:, [int(index), *kept]], rowvar=False
            )[0, 1:]
            if np.max(np.abs(correlations)) < float(threshold):
                kept.append(int(index))
        selected = np.asarray(kept)
    return np.sort(selected)


predictions: dict[str, np.ndarray] = {"current_control": control}
detail_rows: list[dict[str, object]] = []
for objective in cfg["objectives"]:
    objective_predictions = {
        policy["name"]: np.empty(len(y), dtype=float)
        for policy in cfg["pool_policies"]
    }
    for held in np.unique(folds):
        train = folds != held
        valid = ~train
        inner_weights = []
        for inner_held in np.unique(folds[train]):
            inner_train = train & (folds != inner_held)
            inner_weights.append(
                fit_caruana(
                    x[inner_train],
                    y[inner_train],
                    folds[inner_train],
                    objective,
                    np.arange(x.shape[1]),
                )
            )
        inner_weights_array = np.vstack(inner_weights)
        for policy in cfg["pool_policies"]:
            method = f"{objective}::{policy['name']}"
            candidates = adjust_pool(policy, inner_weights_array, x[train])
            weights = fit_caruana(
                x[train], y[train], folds[train], objective, candidates
            )
            objective_predictions[policy["name"]][valid] = x[valid] @ weights
            selected = np.flatnonzero(weights)
            effective_n = 1.0 / float(np.square(weights[selected]).sum())
            detail_rows.append(
                {
                    "method": method,
                    "held_fold": int(held),
                    "pool_size": int(len(candidates)),
                    "selected_count": int(len(selected)),
                    "effective_n": effective_n,
                    "pool_members": json.dumps([names[i] for i in candidates]),
                    "weights": json.dumps(
                        {names[i]: float(weights[i]) for i in selected}, sort_keys=True
                    ),
                    "inner_frequency": json.dumps(
                        {
                            names[i]: int(value)
                            for i, value in enumerate((inner_weights_array > 0).sum(axis=0))
                            if value
                        },
                        sort_keys=True,
                    ),
                }
            )
    for policy in cfg["pool_policies"]:
        method = f"{objective}::{policy['name']}"
        predictions[method] = objective_predictions[policy["name"]]

metric_rows = []
calls_by_method = {}
control_calls = fixed_calls(control, folds)
oof_frames = []
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
    & (pl.col("method") == pool_cfg["jeremy_chemeleon_members"][0]["method"])
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
            "mean_delta": float(delta.mean()),
            "ci025": float(np.quantile(delta, 0.025)),
            "ci975": float(np.quantile(delta, 0.975)),
            "probability_positive": float((delta > 0).mean()),
        }
    )

metrics = pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(detail_rows).write_parquet(out / "outer-details.parquet")
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
pl.concat(oof_frames).write_parquet(out / "oof.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "protocol": cfg["protocol"],
            "initial_pool_size": len(names),
            "members": names,
            "fixed_member_cap": None,
            "candidate_created": False,
            "uploaded": False,
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
