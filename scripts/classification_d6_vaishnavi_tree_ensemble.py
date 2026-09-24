"""Matched CYP2D6 reproduction of the public shallow three-tree ensemble."""

import json
import os
from pathlib import Path

import joblib
import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

from cyp_challenge.tdi_alerts import TDI_ALERTS, alert_matrix

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
    tp = prevalence * tpr
    fn = prevalence * (1 - tpr)
    fp = (1 - prevalence) * fpr
    tn = (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def model_factories(weighting: str, seed: int):
    class_weight = "balanced" if weighting == "balanced" else None
    return {
        "random_forest": RandomForestClassifier(
            n_estimators=600,
            max_depth=8,
            max_features="sqrt",
            class_weight=class_weight,
            n_jobs=3,
            random_state=seed,
        ),
        "extra_trees": ExtraTreesClassifier(
            n_estimators=600,
            max_depth=10,
            max_features="sqrt",
            class_weight=class_weight,
            n_jobs=3,
            random_state=seed,
        ),
        "lightgbm": lgb.LGBMClassifier(
            objective="binary",
            n_estimators=300,
            max_depth=4,
            num_leaves=15,
            learning_rate=0.05,
            class_weight=class_weight,
            n_jobs=3,
            random_state=seed,
            verbosity=-1,
            deterministic=True,
            force_col_wise=True,
        ),
    }


compounds = read("compounds").sort("compound_id")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
labels = read("labels").filter(pl.col("endpoint") == cfg["endpoint"])
assignments = read("assignments")
compound_index = {compound_id: row for row, compound_id in enumerate(compounds["compound_id"])}
test_index = np.asarray([compound_index[value] for value in test["compound_id"]])

molecules = [Chem.MolFromSmiles(value) for value in compounds["standardized_smiles"]]
if any(molecule is None for molecule in molecules):
    raise ValueError("Invalid standardized SMILES")
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan = np.stack([generator.GetFingerprintAsNumPy(molecule) for molecule in molecules]).astype(
    np.float32
)
alert_counts, alert_names = alert_matrix(compounds["standardized_smiles"].to_list())
category_count = len(TDI_ALERTS)
binary_alerts = np.column_stack(
    [
        (alert_counts[:, :category_count] > 0).astype(np.float32),
        alert_counts[:, -2].astype(np.float32),
    ]
)
binary_alert_names = [*alert_names[:category_count], "any_tdi_alert"]
assert binary_alerts.shape[1] == 13

views = {
    "morgan": morgan,
    "morgan_alert13": np.column_stack([morgan, binary_alerts]),
}
pl.DataFrame(
    {"compound_id": compounds["compound_id"]}
    | {f"alert_{name}": binary_alerts[:, index] for index, name in enumerate(binary_alert_names)}
).write_parquet(out / "alert13-features.parquet")

oof_rows = []
test_rows = []
metric_rows = []
model_rows = []
for strategy in cfg["strategies"]:
    split = assignments.filter(pl.col("strategy") == strategy).select(
        "compound_id", "fold", "group"
    )
    data = labels.join(split, on="compound_id", validate="1:1").sort("compound_id")
    data_index = np.asarray([compound_index[value] for value in data["compound_id"]])
    y = data["value"].to_numpy().astype(int)
    folds = data["fold"].to_numpy().astype(int)

    predictions: dict[str, np.ndarray] = {}
    fold_test: dict[str, list[np.ndarray]] = {}
    for view_name, features in views.items():
        for weighting in cfg["weightings"]:
            prefix = f"{view_name}::{weighting}"
            member_names = [
                f"{prefix}::random_forest",
                f"{prefix}::extra_trees",
                f"{prefix}::lightgbm",
            ]
            for name in member_names:
                predictions[name] = np.empty(len(data), dtype=float)
                fold_test[name] = []
            for fold in range(5):
                train = folds != fold
                validation = ~train
                seed = int(cfg["seed"] + fold)
                for family, model in model_factories(weighting, seed).items():
                    name = f"{prefix}::{family}"
                    model.fit(features[data_index[train]], y[train])
                    predictions[name][validation] = model.predict_proba(
                        features[data_index[validation]]
                    )[:, 1]
                    fold_test[name].append(model.predict_proba(features[test_index])[:, 1])
                    joblib.dump(
                        model,
                        out / f"{strategy}-{view_name}-{weighting}-{family}-fold{fold}.joblib",
                        compress=3,
                    )
                    model_rows.append(
                        {
                            "endpoint": cfg["endpoint"],
                            "strategy": strategy,
                            "view": view_name,
                            "weighting": weighting,
                            "family": family,
                            "fold": fold,
                            "seed": seed,
                            "n_train": int(train.sum()),
                            "n_validation": int(validation.sum()),
                            "dimensions": int(features.shape[1]),
                        }
                    )
                print(f"{strategy} {view_name} {weighting} fold {fold}", flush=True)
            ensemble_name = f"{prefix}::equal3"
            predictions[ensemble_name] = np.mean(
                [predictions[name] for name in member_names], axis=0
            )
            fold_test[ensemble_name] = [
                np.mean(
                    [fold_test[name][fold] for name in member_names],
                    axis=0,
                )
                for fold in range(5)
            ]

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
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
write_model_rows = pl.DataFrame(model_rows)
write_model_rows.write_parquet(out / "models.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "endpoint": cfg["endpoint"],
            "source_recipe": {
                "random_forest": "600 trees, max_depth=8",
                "extra_trees": "600 trees, max_depth=10",
                "lightgbm": "300 estimators, max_depth=4, num_leaves=15, learning_rate=0.05",
                "ensemble": "equal probability mean",
            },
            "alert13": binary_alert_names,
            "weightings": cfg["weightings"],
            "selection": "No model or weighting selected before both split results.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics.write_csv(), flush=True)
