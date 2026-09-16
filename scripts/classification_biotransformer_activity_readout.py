"""Evaluate fold-safe product activity summaries for TDI classification."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def fixed_rate_calls(scores: np.ndarray, folds: np.ndarray, rate: float) -> np.ndarray:
    calls = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(rate * len(rows))))
        calls[rows[np.argsort(-scores[rows], kind="stable")[:count]]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray, prevalence: float) -> float:
    positive = y == 1
    tpr, fpr = float(calls[positive].mean()), float(calls[~positive].mean())
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


compounds = read("compounds").sort("compound_id")
labels = read("labels")
assignments = read("assignments")
activity = read("activity_features")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
id_to_index = {compound_id: index for index, compound_id in enumerate(compounds["compound_id"])}
test_indices = np.asarray([id_to_index[compound_id] for compound_id in test["compound_id"]])
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in compounds["standardized_smiles"]
    ]
).astype(np.float32)

metric_rows = []
oof_rows = []
test_rows = []
for endpoint in cfg["endpoints"]:
    enzyme = endpoint.split("_", 1)[0].lower()
    settings = cfg["endpoint_settings"][endpoint]
    endpoint_labels = labels.filter(pl.col("endpoint") == endpoint).select(
        "compound_id", "value"
    )
    for strategy in cfg["strategies"]:
        fold_table = assignments.filter(pl.col("strategy") == strategy).select(
            "compound_id", "fold", "group"
        )
        data = endpoint_labels.join(fold_table, on="compound_id", validate="1:1").sort(
            "compound_id"
        )
        data_indices = np.asarray([id_to_index[compound_id] for compound_id in data["compound_id"]])
        folds = data["fold"].to_numpy().astype(int)
        y = data["value"].to_numpy().astype(int)
        predictions = {method: np.empty(data.height) for method in cfg["methods"]}
        test_predictions = {method: [] for method in cfg["methods"]}
        for fold in range(5):
            fold_features = (
                activity.filter(
                    (pl.col("strategy") == strategy) & (pl.col("model_fold") == fold)
                )
                .drop("strategy", "model_fold")
                .sort("compound_id")
            )
            if not fold_features.select("compound_id").equals(
                compounds.select("compound_id")
            ):
                raise ValueError("Activity feature grid does not match compounds")
            all_columns = [column for column in fold_features.columns if column != "compound_id"]
            endpoint_columns = [column for column in all_columns if column.startswith(enzyme)]
            all_delta_columns = [column for column in all_columns if column.endswith("_delta")]
            endpoint_delta_columns = [
                column for column in endpoint_columns if column.endswith("_delta")
            ]
            blocks = {
                "morgan_lgbm": morgan,
                "activity_endpoint_lgbm": fold_features.select(endpoint_columns).to_numpy(),
                "morgan_activity_endpoint_lgbm": np.column_stack(
                    [morgan, fold_features.select(endpoint_columns).to_numpy()]
                ),
                "morgan_activity_endpoint_delta_lgbm": np.column_stack(
                    [morgan, fold_features.select(endpoint_delta_columns).to_numpy()]
                ),
                "morgan_activity_all_lgbm": np.column_stack(
                    [morgan, fold_features.select(all_columns).to_numpy()]
                ),
                "morgan_activity_all_delta_lgbm": np.column_stack(
                    [morgan, fold_features.select(all_delta_columns).to_numpy()]
                ),
            }
            train_rows = np.flatnonzero(folds != fold)
            heldout_rows = np.flatnonzero(folds == fold)
            for method in cfg["methods"]:
                matrix = blocks[method].astype(np.float32)
                model = lgb.LGBMClassifier(
                    **cfg["lgbm"], verbosity=-1, n_jobs=1, random_state=cfg["seed"] + fold
                )
                model.fit(matrix[data_indices[train_rows]], y[train_rows])
                predictions[method][heldout_rows] = model.predict_proba(
                    matrix[data_indices[heldout_rows]]
                )[:, 1]
                test_predictions[method].append(model.predict_proba(matrix[test_indices])[:, 1])
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
            averaged_test = np.mean(test_predictions[method], axis=0)
            for row, score in enumerate(averaged_test):
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
            "selection": cfg["selection"],
            "methods": cfg["methods"],
            "activity_feature_count": len(activity.columns) - 3,
            "outer_safe_routing": True,
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
