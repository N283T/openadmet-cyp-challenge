"""Dual-split controlled readouts for atom-derived SMARTCyp features."""

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


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


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
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def impute_scale(
    train: np.ndarray, query: np.ndarray, scale: bool
) -> tuple[np.ndarray, np.ndarray]:
    medians = np.nanmedian(train, axis=0)
    train = np.where(np.isfinite(train), train, medians)
    query = np.where(np.isfinite(query), query, medians)
    if not scale:
        return train, query
    scaler = StandardScaler().fit(train)
    return scaler.transform(train), scaler.transform(query)


def predict(method: str, x_train: np.ndarray, y: np.ndarray, x_query: np.ndarray) -> np.ndarray:
    if method.endswith("_lr"):
        model = LogisticRegression(
            C=cfg["lr_c"], solver="liblinear", max_iter=4000, random_state=cfg["seed"]
        ).fit(x_train, y)
    else:
        model = lgb.LGBMClassifier(
            n_estimators=cfg["lgbm"]["n_estimators"],
            learning_rate=cfg["lgbm"]["learning_rate"],
            num_leaves=cfg["lgbm"]["num_leaves"],
            min_child_samples=cfg["lgbm"]["min_child_samples"],
            reg_lambda=cfg["lgbm"]["reg_lambda"],
            verbosity=-1,
            n_jobs=1,
            random_state=cfg["seed"],
        ).fit(x_train, y)
    return model.predict_proba(x_query)[:, 1]


compounds = read("compounds").sort("compound_id")
labels = read("labels")
assignments = read("assignments")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
all_ids = compounds.select("compound_id")
id_to_index = {value: index for index, value in enumerate(all_ids["compound_id"])}
test_index = np.asarray([id_to_index[value] for value in test["compound_id"]])
molecules = [Chem.MolFromSmiles(value) for value in compounds["standardized_smiles"]]
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan = np.stack([generator.GetFingerprintAsNumPy(molecule) for molecule in molecules]).astype(
    np.float32
)
smartcyp_table = all_ids.join(read("smartcyp"), on="compound_id", validate="1:1")
all_smartcyp_columns = [column for column in smartcyp_table.columns if column != "compound_id"]
core_smartcyp_columns = [
    column for column in all_smartcyp_columns if "tdi_alert" not in column
]
smartcyp_all = smartcyp_table.select(all_smartcyp_columns).to_numpy().astype(np.float32)
smartcyp_core = smartcyp_table.select(core_smartcyp_columns).to_numpy().astype(np.float32)
feature_sets = {
    "morgan": morgan,
    "smartcyp_core": smartcyp_core,
    "smartcyp_all": smartcyp_all,
    "morgan_smartcyp_core": np.column_stack([morgan, smartcyp_core]),
    "morgan_smartcyp_all": np.column_stack([morgan, smartcyp_all]),
}

oof_rows = []
test_rows = []
metric_rows = []
for endpoint in cfg["endpoints"]:
    settings = cfg["endpoint_settings"][endpoint]
    endpoint_labels = labels.filter(pl.col("endpoint") == endpoint).select("compound_id", "value")
    for strategy in cfg["strategies"]:
        fold_table = assignments.filter(pl.col("strategy") == strategy).select(
            "compound_id", "fold", "group"
        )
        data = endpoint_labels.join(fold_table, on="compound_id", validate="1:1").sort(
            "compound_id"
        )
        indices = np.asarray([id_to_index[value] for value in data["compound_id"]])
        y = data["value"].to_numpy().astype(int)
        folds = data["fold"].to_numpy().astype(int)
        predictions = {method: np.empty(len(data)) for method in cfg["methods"]}
        test_predictions = {method: [] for method in cfg["methods"]}
        for fold in range(5):
            train_rows = np.flatnonzero(folds != fold)
            heldout_rows = np.flatnonzero(folds == fold)
            train_indices = indices[train_rows]
            heldout_indices = indices[heldout_rows]
            for method in cfg["methods"]:
                feature_name = method.removesuffix("_lr").removesuffix("_lgbm")
                matrix = feature_sets[feature_name]
                x_train, x_heldout = impute_scale(
                    matrix[train_indices], matrix[heldout_indices], method.endswith("_lr")
                )
                _, x_test = impute_scale(
                    matrix[train_indices], matrix[test_index], method.endswith("_lr")
                )
                predictions[method][heldout_rows] = predict(
                    method, x_train, y[train_rows], x_heldout
                )
                test_predictions[method].append(predict(method, x_train, y[train_rows], x_test))
            print(f"{endpoint} {strategy} fold {fold} complete", flush=True)

        for method, probability in predictions.items():
            calls = fixed_rate_calls(probability, folds, settings["call_rate"])
            metric_rows.append(
                {
                    "endpoint": endpoint,
                    "strategy": strategy,
                    "method": method,
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
                        "method": method,
                        "prediction": float(score),
                        "call": bool(calls[row]),
                    }
                )
            for row, score in enumerate(np.mean(test_predictions[method], axis=0)):
                test_rows.append(
                    {
                        "compound_id": test["compound_id"][row],
                        "endpoint": endpoint,
                        "strategy": strategy,
                        "fold": -1,
                        "method": method,
                        "prediction": float(score),
                    }
                )

metrics = pl.DataFrame(metric_rows).sort("endpoint", "strategy", "method")
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(oof_rows).write_parquet(out / "oof.parquet")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "smartcyp_columns": all_smartcyp_columns,
            "core_columns": core_smartcyp_columns,
            "methods": cfg["methods"],
            "strategies": cfg["strategies"],
            "selection": "Prespecified controlled comparison; no submission selection.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
