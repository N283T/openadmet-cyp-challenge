"""Screen TabPFN 3.5-fast on selected CYP2D6 classification feature pools."""

import importlib.metadata
import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
import torch
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from tabpfn import TabPFNClassifier
from tabpfn.constants import ModelVersion
from threadpoolctl import threadpool_limits

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
torch.set_num_threads(cfg["num_threads"])
threadpool_limits(cfg["num_threads"])


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def load_block(name: str) -> tuple[list[str], np.ndarray]:
    spec = cfg["blocks"][name]
    table = pl.read_parquet(inp / spec["file"]).sort("compound_id")
    if table["compound_id"].to_list() != canonical_ids:
        raise ValueError(f"Compound order mismatch in block {name}")
    columns = [
        column
        for column in table.columns
        if column != "compound_id" and any(column.startswith(prefix) for prefix in spec["prefixes"])
    ]
    values = table.select(columns).to_numpy().astype(np.float32)
    values[~np.isfinite(values)] = np.nan
    return [f"{name}::{column}" for column in columns], values


def fixed_rate_calls(scores: np.ndarray, folds: np.ndarray) -> np.ndarray:
    calls = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(positions))))
        calls[positions[np.argsort(-scores[positions], kind="stable")[:count]]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray) -> float:
    positive = y == 1
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    prevalence = cfg["target_prior"]
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def finite_matrix(values: np.ndarray, training: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    medians = np.nanmedian(values[training], axis=0)
    medians[~np.isfinite(medians)] = 0.0
    result = np.where(np.isfinite(values), values, medians).astype(np.float32)
    if not np.isfinite(result).all():
        raise ValueError("Non-finite values remain after training-median imputation")
    return result, medians


def select_features(
    values: np.ndarray,
    y: np.ndarray,
    training: np.ndarray,
    names: list[str],
    fold: int,
    top_k: int,
) -> tuple[np.ndarray, np.ndarray]:
    selector = lgb.LGBMClassifier(
        **cfg["selector"],
        objective="binary",
        class_weight="balanced",
        n_jobs=cfg["num_threads"],
        random_state=cfg["seed"] + fold,
        verbosity=-1,
        deterministic=True,
        force_col_wise=True,
    )
    x, _ = finite_matrix(values, training)
    selector.fit(x[training], y[training])
    gain = selector.booster_.feature_importance(importance_type="gain")
    order = np.lexsort((np.asarray(names), -gain))
    return order[: min(top_k, len(order))], gain


def make_model(fold: int) -> TabPFNClassifier:
    return TabPFNClassifier.create_default_for_version(
        ModelVersion(cfg["tabpfn_version"]),
        n_estimators=cfg["tabpfn_estimators"],
        model_path=str(inp / cfg["tabpfn_weights"]),
        device="cuda",
        random_state=cfg["seed"] + fold,
        auto_scale_n_estimators=False,
        n_preprocessing_jobs=1,
        memory_saving_mode="auto",
        show_progress_bar=False,
        ignore_pretraining_limits=True,
        fit_mode="fit_preprocessors",
    )


identity = (
    pl.read_parquet(inp / cfg["blocks"]["morgan"]["file"]).select("compound_id").sort("compound_id")
)
canonical_ids = identity["compound_id"].to_list()
id_to_index = {compound_id: row for row, compound_id in enumerate(canonical_ids)}
blocks = {name: load_block(name) for name in cfg["blocks"]}
labels = read("labels").filter(pl.col("endpoint") == cfg["endpoint"])
assignments = read("assignments")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
test_index = np.asarray([id_to_index[value] for value in test["compound_id"]])

candidates = {}
candidate_names = {}
for candidate in cfg["candidates"]:
    names = []
    matrices = []
    for block in candidate["blocks"]:
        block_names, block_values = blocks[block]
        names.extend(block_names)
        matrices.append(block_values)
    candidate_names[candidate["name"]] = names
    candidates[candidate["name"]] = np.column_stack(matrices)

oof_rows = []
test_rows = []
metric_rows = []
fit_rows = []
selection_rows = []
all_pool_rows = []
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
        name = candidate["name"]
        values = candidates[name]
        names = candidate_names[name]
        predictions[name] = np.empty(len(data), dtype=float)
        fold_test[name] = []
        for fold in range(5):
            training = folds != fold
            validation = ~training
            selected, gain = select_features(
                values[data_index],
                y,
                training,
                names,
                fold,
                candidate["importance_top_k"],
            )
            selected_values, medians = finite_matrix(values[:, selected], data_index[training])
            model = make_model(fold)
            model.fit(selected_values[data_index[training]], y[training])
            pool_valid = pool_folds == fold
            pool_scores = model.predict_proba(selected_values[pool_index[pool_valid]])[:, 1]
            if not np.isfinite(pool_scores).all():
                raise ValueError(f"Nonfinite all-pool score: {strategy} {name} {fold}")
            pool_predictions = (
                pool.filter(pl.Series(pool_valid))
                .select("compound_id")
                .with_columns(pl.Series("score", pool_scores))
            )
            labelled_predictions = (
                data.filter(pl.Series(validation))
                .select("compound_id")
                .join(pool_predictions, on="compound_id", validate="1:1", maintain_order="left")
            )
            if labelled_predictions["score"].null_count():
                raise ValueError(
                    f"Labelled prediction missing from all-pool batch: {strategy} {fold}"
                )
            predictions[name][validation] = labelled_predictions["score"].to_numpy()
            all_pool_rows.append(
                pool.filter(pl.Series(pool_valid)).with_columns(
                    pl.lit(cfg["endpoint"]).alias("endpoint"),
                    pl.lit(strategy).alias("split_id"),
                    pl.lit(name).alias("model"),
                    pl.Series("score", pool_scores),
                )
            )
            fold_test[name].append(model.predict_proba(selected_values[test_index])[:, 1])
            fit_rows.append(
                {
                    "strategy": strategy,
                    "candidate": name,
                    "fold": fold,
                    "input_dimensions": values.shape[1],
                    "selected_dimensions": len(selected),
                    "training": int(training.sum()),
                    "validation": int(validation.sum()),
                }
            )
            for rank, feature_index in enumerate(selected, start=1):
                selection_rows.append(
                    {
                        "strategy": strategy,
                        "candidate": name,
                        "fold": fold,
                        "rank": rank,
                        "feature": names[feature_index],
                        "selector_gain": float(gain[feature_index]),
                        "training_median": float(medians[rank - 1]),
                    }
                )
            print(json.dumps(fit_rows[-1]), flush=True)
            del model
            torch.cuda.empty_cache()

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
        calls = fixed_rate_calls(probability, folds)
        metric_rows.append(
            {
                "strategy": strategy,
                "method": method,
                "AUROC": float(roc_auc_score(y, probability)),
                "AP": float(average_precision_score(y, probability)),
                "MCC_observed": float(matthews_corrcoef(y, calls)),
                "MCC_target_prior": expected_mcc(y, calls),
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
            averaged = np.mean(fold_test[method], axis=0)
            count = int(round(cfg["call_rate"] * len(test)))
            test_calls = np.zeros(len(test), dtype=bool)
            test_calls[np.argsort(-averaged, kind="stable")[:count]] = True
            for row, prediction in enumerate(averaged):
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
pl.DataFrame(fit_rows).write_csv(out / "fits.csv")
pl.DataFrame(selection_rows).write_parquet(out / "feature-selections.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "endpoint": cfg["endpoint"],
            "tabpfn_package": importlib.metadata.version("tabpfn"),
            "tabpfn_model_version": cfg["tabpfn_version"],
            "tabpfn_estimators": cfg["tabpfn_estimators"],
            "selection": cfg["selection"],
            "candidate_created": False,
            "uploaded": False,
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
