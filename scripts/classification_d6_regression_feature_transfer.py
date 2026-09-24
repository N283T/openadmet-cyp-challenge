"""Transfer regression feature pools into fold-safe CYP2D6 classifiers."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def load_block(name: str) -> tuple[list[str], dict[int | None, np.ndarray]]:
    spec = cfg["blocks"][name]
    table = pl.read_parquet(inp / spec["file"])
    columns = [
        column
        for column in table.columns
        if column != "compound_id" and any(column.startswith(prefix) for prefix in spec["prefixes"])
    ]
    if not columns:
        raise ValueError(f"No selected columns in block {name}")
    fold_column = spec.get("outer_fold_column")
    folds = range(5) if fold_column else [None]
    matrices = {}
    for fold in folds:
        frame = table.filter(pl.col(fold_column) == fold) if fold_column else table
        frame = frame.sort("compound_id")
        if frame["compound_id"].to_list() != canonical_ids:
            raise ValueError(f"Compound order mismatch in block {name}, fold {fold}")
        matrices[fold] = frame.select(columns).to_numpy().astype(np.float32)
    return [f"{name}::{column}" for column in columns], matrices


def fixed_rate_calls(scores: np.ndarray, folds: np.ndarray, rate: float) -> np.ndarray:
    calls = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = max(1, int(round(rate * len(positions))))
        order = np.argsort(-scores[positions], kind="stable")
        calls[positions[order[:count]]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray, prevalence: float) -> float:
    positive = y == 1
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    tp = prevalence * tpr
    fn = prevalence * (1 - tpr)
    fp = (1 - prevalence) * fpr
    tn = (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def make_model(weighting: str, seed: int, n_estimators: int | None = None):
    params = cfg["lightgbm"] | {
        "objective": "binary",
        "class_weight": "balanced" if weighting == "balanced" else None,
        "n_jobs": 4,
        "random_state": seed,
        "verbosity": -1,
        "deterministic": True,
        "force_col_wise": True,
    }
    if n_estimators is not None:
        params["n_estimators"] = n_estimators
    return lgb.LGBMClassifier(**params)


identity_file = next(iter(cfg["blocks"].values()))["file"]
identity = pl.read_parquet(inp / identity_file).select("compound_id").sort("compound_id")
canonical_ids = identity["compound_id"].to_list()
id_to_index = {compound_id: row for row, compound_id in enumerate(canonical_ids)}
blocks = {name: load_block(name) for name in cfg["blocks"]}
labels = read("labels").filter(pl.col("endpoint") == cfg["endpoint"])
assignments = read("assignments")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
test_index = np.asarray([id_to_index[value] for value in test["compound_id"]])

candidate_names = {}
for candidate in cfg["candidates"]:
    names = []
    for block in candidate["blocks"]:
        block_names, _ = blocks[block]
        names.extend(block_names)
    candidate_names[candidate["name"]] = names


def candidate_matrix(candidate: dict, fold: int) -> np.ndarray:
    matrices = []
    for block in candidate["blocks"]:
        _, block_values = blocks[block]
        matrices.append(block_values.get(fold, block_values.get(None)))
    return np.column_stack(matrices)


oof_rows = []
all_pool_rows = []
test_rows = []
metric_rows = []
selection_rows = []
for strategy in cfg["strategies"]:
    split = assignments.filter(pl.col("strategy") == strategy).select(
        "compound_id", "fold", "group"
    )
    if split.height != 6146 or split["compound_id"].n_unique() != 6146:
        raise ValueError(f"Incomplete molecular pool for {strategy}")
    pool = split.join(
        labels.select("compound_id", pl.col("value").cast(pl.Int64).alias("y_true")),
        on="compound_id",
        how="left",
        validate="1:1",
    ).sort("compound_id")
    pool_index = np.asarray([id_to_index[value] for value in pool["compound_id"]])
    pool_folds = pool["fold"].to_numpy().astype(int)
    data = labels.join(split, on="compound_id", validate="1:1").sort("compound_id")
    data_index = np.asarray([id_to_index[value] for value in data["compound_id"]])
    y = data["value"].to_numpy().astype(int)
    folds = data["fold"].to_numpy().astype(int)
    predictions = {}
    fold_test = {}
    for candidate in cfg["candidates"]:
        candidate_name = candidate["name"]
        feature_names = candidate_names[candidate_name]
        for weighting in cfg["weightings"]:
            method = f"{candidate_name}::{weighting}"
            predictions[method] = np.empty(len(data), dtype=float)
            fold_test[method] = []
            for fold in range(5):
                features = candidate_matrix(candidate, fold)
                train = folds != fold
                validation = ~train
                selected = np.arange(features.shape[1])
                if candidate.get("importance_top_k"):
                    selector = make_model(weighting, cfg["seed"] + fold, n_estimators=200)
                    selector.fit(features[data_index[train]], y[train])
                    importance = selector.feature_importances_
                    selected = np.argsort(-importance, kind="stable")[
                        : candidate["importance_top_k"]
                    ]
                model = make_model(weighting, cfg["seed"] + fold)
                model.fit(features[data_index[train]][:, selected], y[train])
                predictions[method][validation] = model.predict_proba(
                    features[data_index[validation]][:, selected]
                )[:, 1]
                pool_validation = pool_folds == fold
                pool_scores = model.predict_proba(
                    features[pool_index[pool_validation]][:, selected]
                )[:, 1]
                if not np.isfinite(pool_scores).all():
                    raise ValueError(f"Nonfinite all-pool prediction: {strategy} {method} {fold}")
                all_pool_rows.append(
                    pool.filter(pl.Series(pool_validation)).with_columns(
                        pl.lit(cfg["endpoint"]).alias("endpoint"),
                        pl.lit(strategy).alias("split_id"),
                        pl.lit(method).alias("model"),
                        pl.Series("score", pool_scores),
                    )
                )
                fold_test[method].append(
                    model.predict_proba(features[test_index][:, selected])[:, 1]
                )
                selected_names = [feature_names[index] for index in selected]
                selection_rows.append(
                    {
                        "endpoint": cfg["endpoint"],
                        "strategy": strategy,
                        "candidate": candidate_name,
                        "weighting": weighting,
                        "fold": fold,
                        "input_dimensions": int(features.shape[1]),
                        "selected_dimensions": int(len(selected)),
                        "selected_features": json.dumps(selected_names),
                    }
                )
                print(f"{strategy} {method} fold {fold}", flush=True)

    if strategy == "umap":
        current = (
            read("current_bank")
            .filter(
                (pl.col("endpoint") == cfg["endpoint"]) & (pl.col("model") == "current_fixed_ratio")
            )
            .select("compound_id", pl.col("prediction").alias("current_prediction"))
        )
        control = data.select("compound_id").join(
            current, on="compound_id", validate="1:1", maintain_order="left"
        )
        predictions["current_fixed_ratio"] = control["current_prediction"].to_numpy()

    for method, probability in predictions.items():
        calls = fixed_rate_calls(probability, folds, cfg["call_rate"])
        metric_rows.append(
            {
                "endpoint": cfg["endpoint"],
                "strategy": strategy,
                "method": method,
                "AUROC": float(roc_auc_score(y, probability)),
                "AP": float(average_precision_score(y, probability)),
                "MCC_observed": float(matthews_corrcoef(y, calls)),
                "MCC_target_prior": expected_mcc(y, calls, cfg["target_prior"]),
                "positive_calls": int(calls.sum()),
            }
        )
        for row, prediction in enumerate(probability):
            oof_rows.append(
                {
                    "compound_id": data["compound_id"][row],
                    "endpoint": cfg["endpoint"],
                    "split_id": strategy,
                    "fold": int(folds[row]),
                    "group": data["group"][row],
                    "y_true": int(y[row]),
                    "method": method,
                    "prediction": float(prediction),
                    "call": bool(calls[row]),
                }
            )
        if method in fold_test:
            averaged_test = np.mean(fold_test[method], axis=0)
            test_calls = np.zeros(len(test), dtype=bool)
            count = int(round(cfg["call_rate"] * len(test)))
            test_calls[np.argsort(-averaged_test, kind="stable")[:count]] = True
            for row, prediction in enumerate(averaged_test):
                test_rows.append(
                    {
                        "compound_id": test["compound_id"][row],
                        "endpoint": cfg["endpoint"],
                        "split_id": strategy,
                        "fold": -1,
                        "method": method,
                        "prediction": float(prediction),
                        "call": bool(test_calls[row]),
                    }
                )

metrics = pl.DataFrame(metric_rows).sort("strategy", "MCC_target_prior", descending=[False, True])
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(oof_rows).write_parquet(out / "oof.parquet")
all_pool = pl.concat(all_pool_rows).sort("split_id", "model", "compound_id")
for (split_id, model), frame in all_pool.group_by("split_id", "model"):
    if frame.height != 6146 or frame["compound_id"].n_unique() != 6146:
        raise ValueError(f"All-pool coverage failed: {split_id} {model}")
all_pool.write_parquet(out / "all-pool-oof.parquet")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
pl.DataFrame(selection_rows).write_parquet(out / "feature-selections.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "endpoint": cfg["endpoint"],
            "candidate_dimensions": {
                candidate["name"]: int(candidate_matrix(candidate, 0).shape[1])
                for candidate in cfg["candidates"]
            },
            "selection": "Importance selection is repeated inside each outer-training fold.",
            "test_policy": "49 highest mean fold probabilities; no threshold tuning.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics.write_csv(), flush=True)
