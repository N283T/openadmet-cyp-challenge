"""Build heterogeneous readouts around the anchored Jeremy/CheMeleon backbone."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from catboost import CatBoostClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from xgboost import XGBClassifier

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def fixed_calls(scores: np.ndarray, folds: np.ndarray) -> np.ndarray:
    calls = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(rows))))
        calls[rows[np.argsort(-scores[rows], kind="stable")[:count]]] = True
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


identity_file = next(iter(cfg["blocks"].values()))["file"]
identity = pl.read_parquet(inp / identity_file).select("compound_id").sort("compound_id")
canonical_ids = identity["compound_id"].to_list()
id_to_index = {value: index for index, value in enumerate(canonical_ids)}


def load_block(name: str) -> tuple[list[str], np.ndarray]:
    spec = cfg["blocks"][name]
    table = pl.read_parquet(inp / spec["file"]).sort("compound_id")
    if table["compound_id"].to_list() != canonical_ids:
        raise ValueError(f"Compound order mismatch in {name}")
    columns = [
        column
        for column in table.columns
        if column != "compound_id" and any(column.startswith(prefix) for prefix in spec["prefixes"])
    ]
    if not columns:
        raise ValueError(f"No features selected for {name}")
    values = table.select(columns).to_numpy().astype(np.float32)
    values[~np.isfinite(values)] = np.nan
    return [f"{name}::{column}" for column in columns], values


blocks = {name: load_block(name) for name in cfg["blocks"]}
candidate_data = {}
for candidate in cfg["candidates"]:
    names, matrices = [], []
    for block in candidate["blocks"]:
        block_names, block_values = blocks[block]
        names.extend(block_names)
        matrices.append(block_values)
    candidate_data[candidate["name"]] = (names, np.column_stack(matrices))

labels = read("labels").filter(pl.col("endpoint") == cfg["endpoint"])
assignments = read("assignments")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
test_index = np.asarray([id_to_index[value] for value in test["compound_id"]])


def selector(seed: int):
    return lgb.LGBMClassifier(
        **cfg["selector"],
        objective="binary",
        class_weight="balanced",
        n_jobs=cfg["num_threads"],
        random_state=seed,
        verbosity=-1,
        deterministic=True,
        force_col_wise=True,
    )


def make_readout(name: str, seed: int):
    if name == "catboost":
        return CatBoostClassifier(
            iterations=220,
            depth=5,
            learning_rate=0.035,
            l2_leaf_reg=5.0,
            loss_function="Logloss",
            auto_class_weights="Balanced",
            random_seed=seed,
            thread_count=cfg["num_threads"],
            task_type="GPU",
            devices="0",
            verbose=False,
            allow_writing_files=False,
        )
    if name == "xgboost":
        return XGBClassifier(
            n_estimators=240,
            max_depth=4,
            learning_rate=0.03,
            min_child_weight=4,
            subsample=0.85,
            colsample_bytree=0.7,
            reg_lambda=5.0,
            objective="binary:logistic",
            eval_metric="logloss",
            n_jobs=cfg["num_threads"],
            random_state=seed,
            device="cuda",
        )
    if name == "extra_trees":
        return ExtraTreesClassifier(
            n_estimators=400,
            min_samples_leaf=2,
            max_features="sqrt",
            class_weight="balanced",
            n_jobs=cfg["num_threads"],
            random_state=seed,
        )
    raise ValueError(name)


oof_rows, test_rows, metric_rows, selection_rows, all_pool_rows = [], [], [], [], []
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
    for candidate in cfg["candidates"]:
        candidate_name = candidate["name"]
        feature_names, values = candidate_data[candidate_name]
        predictions = {readout: np.empty(len(data), dtype=float) for readout in cfg["readouts"]}
        fold_test = {readout: [] for readout in cfg["readouts"]}
        for fold in range(5):
            train = folds != fold
            valid = ~train
            train_rows = data_index[train]
            finite = np.where(np.isfinite(values), values, np.nan)
            medians = np.nanmedian(finite[train_rows], axis=0)
            medians[~np.isfinite(medians)] = 0.0
            filled = np.where(np.isfinite(finite), finite, medians).astype(np.float32)
            choose = selector(cfg["seed"] + fold)
            choose.fit(filled[train_rows], y[train])
            gain = choose.booster_.feature_importance(importance_type="gain")
            order = np.lexsort((np.asarray(feature_names), -gain))
            selected = order[: min(cfg["importance_top_k"], len(order))]
            x_train = filled[train_rows][:, selected]
            x_valid = filled[data_index[valid]][:, selected]
            x_test = filled[test_index][:, selected]
            pool_valid = pool_folds == fold
            x_pool = filled[pool_index[pool_valid]][:, selected]
            offsets = {"catboost": 0, "xgboost": 1, "extra_trees": 2}
            for readout in cfg["readouts"]:
                model = make_readout(readout, cfg["seed"] + 100 * offsets[readout] + fold)
                model.fit(x_train, y[train])
                predictions[readout][valid] = model.predict_proba(x_valid)[:, 1]
                pool_score = model.predict_proba(x_pool)[:, 1]
                if not np.isfinite(pool_score).all():
                    raise ValueError(
                        f"Nonfinite all-pool score: {strategy} {candidate_name} {readout}"
                    )
                all_pool_rows.append(
                    pool.filter(pl.Series(pool_valid)).with_columns(
                        pl.lit(cfg["endpoint"]).alias("endpoint"),
                        pl.lit(strategy).alias("split_id"),
                        pl.lit(f"{candidate_name}::{readout}").alias("model"),
                        pl.Series("score", pool_score),
                    )
                )
                fold_test[readout].append(model.predict_proba(x_test)[:, 1])
            selection_rows.append(
                {
                    "strategy": strategy,
                    "candidate": candidate_name,
                    "fold": fold,
                    "input_dimensions": values.shape[1],
                    "selected_dimensions": len(selected),
                    "selected_features": json.dumps([feature_names[index] for index in selected]),
                }
            )
            print(f"{strategy} {candidate_name} fold {fold}", flush=True)
        for readout, probability in predictions.items():
            method = f"{candidate_name}::{readout}"
            calls = fixed_calls(probability, folds)
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
            averaged = np.mean(fold_test[readout], axis=0)
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
pl.DataFrame(selection_rows).write_parquet(out / "feature-selections.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "endpoint": cfg["endpoint"],
            "protocol": cfg["protocol"],
            "candidate_created": False,
            "uploaded": False,
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
