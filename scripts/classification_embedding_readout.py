"""Controlled dual-split readout for frozen molecular feature blocks."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from sklearn.preprocessing import StandardScaler

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(key):
    return pl.read_parquet(inp / cfg[key])


def fixed_rate_calls(scores, folds, rate):
    calls = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = max(1, int(round(rate * len(positions))))
        calls[positions[np.argsort(-scores[positions], kind="stable")[:count]]] = True
    return calls


def expected_mcc(y, calls, prevalence):
    positive = y == 1
    tpr, fpr = float(calls[positive].mean()), float(calls[~positive].mean())
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def preprocess(train, *queries, scale):
    medians = np.nanmedian(train, axis=0)
    train = np.where(np.isfinite(train), train, medians)
    queries = [np.where(np.isfinite(query), query, medians) for query in queries]
    if not scale:
        return (train, *queries)
    scaler = StandardScaler().fit(train)
    return (scaler.transform(train), *[scaler.transform(query) for query in queries])


def fit_predict(model_name, x_train, y_train, *queries):
    if model_name == "lr":
        model = LogisticRegression(
            C=cfg["lr_c"], solver="liblinear", max_iter=4000, random_state=cfg["seed"]
        )
    else:
        model = lgb.LGBMClassifier(**cfg["lgbm"], verbosity=-1, n_jobs=1, random_state=cfg["seed"])
    model.fit(x_train, y_train)
    return [model.predict_proba(query)[:, 1] for query in queries]


compounds = read("compounds").sort("compound_id")
all_ids = compounds.select("compound_id")
id_to_index = {value: index for index, value in enumerate(all_ids["compound_id"])}
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
test_indices = np.asarray([id_to_index[value] for value in test["compound_id"]])
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
feature_blocks = {
    "morgan": np.stack(
        [
            generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
            for smiles in compounds["standardized_smiles"]
        ]
    ).astype(np.float32)
}
feature_columns = {}
for block, source_spec in cfg["feature_sources"].items():
    filename = source_spec if isinstance(source_spec, str) else source_spec["file"]
    table = all_ids.join(
        pl.read_parquet(inp / filename), on="compound_id", how="left", validate="1:1"
    )
    columns = [column for column in table.columns if column != "compound_id"]
    if not isinstance(source_spec, str) and "columns" in source_spec:
        requested = source_spec["columns"]
        missing = sorted(set(requested) - set(columns))
        if missing:
            raise ValueError(f"Missing columns for feature block {block}: {missing}")
        columns = requested
    if not isinstance(source_spec, str) and "exclude_substrings" in source_spec:
        columns = [
            column
            for column in columns
            if not any(token in column for token in source_spec["exclude_substrings"])
        ]
    feature_blocks[block] = table.select(columns).to_numpy().astype(np.float32)
    feature_columns[block] = columns

matrices = {
    method["name"]: np.column_stack([feature_blocks[block] for block in method["blocks"]])
    for method in cfg["methods"]
}
method_by_name = {method["name"]: method for method in cfg["methods"]}
labels = read("labels")
assignments = read("assignments")
oof_rows, test_rows, metric_rows, all_pool_rows = [], [], [], []
for endpoint in cfg["endpoints"]:
    settings = cfg["endpoint_settings"][endpoint]
    endpoint_labels = labels.filter(pl.col("endpoint") == endpoint).select("compound_id", "value")
    for strategy in cfg["strategies"]:
        fold_table = assignments.filter(pl.col("strategy") == strategy).select(
            "compound_id", "fold", "group"
        )
        if fold_table.height != 6146 or fold_table["compound_id"].n_unique() != 6146:
            raise ValueError(f"Incomplete molecular pool for {strategy}")
        pool = fold_table.join(
            endpoint_labels.select("compound_id", pl.col("value").cast(pl.Int64).alias("y_true")),
            on="compound_id",
            how="left",
            validate="1:1",
        ).sort("compound_id")
        pool_indices = np.asarray([id_to_index[value] for value in pool["compound_id"]])
        pool_folds = pool["fold"].to_numpy().astype(int)
        data = endpoint_labels.join(fold_table, on="compound_id", validate="1:1").sort(
            "compound_id"
        )
        indices = np.asarray([id_to_index[value] for value in data["compound_id"]])
        y = data["value"].to_numpy().astype(int)
        folds = data["fold"].to_numpy().astype(int)
        predictions = {name: np.empty(len(data)) for name in matrices}
        test_predictions = {name: [] for name in matrices}
        for fold in range(5):
            train_rows, heldout_rows = np.flatnonzero(folds != fold), np.flatnonzero(folds == fold)
            for name, matrix in matrices.items():
                method = method_by_name[name]
                pool_valid = pool_folds == fold
                x_train, x_heldout, x_test, x_pool = preprocess(
                    matrix[indices[train_rows]],
                    matrix[indices[heldout_rows]],
                    matrix[test_indices],
                    matrix[pool_indices[pool_valid]],
                    scale=method["model"] == "lr",
                )
                heldout_pred, test_pred, pool_pred = fit_predict(
                    method["model"], x_train, y[train_rows], x_heldout, x_test, x_pool
                )
                predictions[name][heldout_rows] = heldout_pred
                test_predictions[name].append(test_pred)
                if not np.isfinite(pool_pred).all():
                    raise ValueError(f"Nonfinite all-pool score: {endpoint} {strategy} {name}")
                all_pool_rows.append(
                    pool.filter(pl.Series(pool_valid)).with_columns(
                        pl.lit(endpoint).alias("endpoint"),
                        pl.lit(strategy).alias("split_id"),
                        pl.lit(name).alias("model"),
                        pl.Series("score", pool_pred),
                    )
                )
            print(f"{endpoint} {strategy} fold {fold} complete", flush=True)
        for name, probability in predictions.items():
            calls = fixed_rate_calls(probability, folds, settings["call_rate"])
            metric_rows.append(
                {
                    "endpoint": endpoint,
                    "strategy": strategy,
                    "method": name,
                    "AUROC": float(roc_auc_score(y, probability)),
                    "AP": float(average_precision_score(y, probability)),
                    "MCC_observed": float(matthews_corrcoef(y, calls)),
                    "MCC_target_prior": expected_mcc(y, calls, settings["target_prior"]),
                    "calls": int(calls.sum()),
                }
            )
            for row, score in enumerate(probability):
                oof_rows.append(
                    {
                        "compound_id": data["compound_id"][row],
                        "endpoint": endpoint,
                        "strategy": strategy,
                        "fold": int(folds[row]),
                        "group": data["group"][row],
                        "y_true": int(y[row]),
                        "method": name,
                        "prediction": float(score),
                        "call": bool(calls[row]),
                    }
                )
            for row, score in enumerate(np.mean(test_predictions[name], axis=0)):
                test_rows.append(
                    {
                        "compound_id": test["compound_id"][row],
                        "endpoint": endpoint,
                        "strategy": strategy,
                        "fold": -1,
                        "method": name,
                        "prediction": float(score),
                    }
                )

metrics = pl.DataFrame(metric_rows).sort("endpoint", "strategy", "method")
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(oof_rows).write_parquet(out / "oof.parquet")
all_pool = pl.concat(all_pool_rows).sort("endpoint", "split_id", "model", "compound_id")
for (endpoint, split_id, model), frame in all_pool.group_by("endpoint", "split_id", "model"):
    if frame.height != 6146 or frame["compound_id"].n_unique() != 6146:
        raise ValueError(f"All-pool coverage failed: {endpoint} {split_id} {model}")
all_pool.write_parquet(out / "all-pool-oof.parquet")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "feature_columns": feature_columns,
            "methods": cfg["methods"],
            "strategies": cfg["strategies"],
            "selection": cfg["selection"],
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
