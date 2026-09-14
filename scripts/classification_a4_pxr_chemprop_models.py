"""Compare model families on the fixed strict-PXR-Chemprop embedding."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC, LinearSVC

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def fixed_calls(scores: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(positions))))
        result[positions[np.argsort(-scores[positions], kind="stable")[:count]]] = True
    return result


def expected_mcc(y: np.ndarray, selected: np.ndarray) -> float:
    positive = y == 1
    tpr = float(selected[positive].mean())
    fpr = float(selected[~positive].mean())
    prevalence = cfg["target_prior"]
    tp = prevalence * tpr
    fn = prevalence * (1 - tpr)
    fp = (1 - prevalence) * fpr
    tn = (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def model_specs(seed: int):
    result = []
    for c_value in [0.01, 0.05, 0.2, 1.0]:
        result.append(
            (
                f"logistic_c{c_value:g}",
                LogisticRegression(C=c_value, solver="liblinear", max_iter=4000, random_state=seed),
                "scaled",
            )
        )
    for c_value in [0.001, 0.01, 0.1, 1.0]:
        result.append((f"linear_svm_c{c_value:g}", LinearSVC(C=c_value), "scaled"))
    for c_value in [0.1, 1.0, 10.0]:
        result.append((f"rbf_svm_c{c_value:g}", SVC(C=c_value), "scaled"))
        result.append((f"pca64_rbf_svm_c{c_value:g}", SVC(C=c_value), "pca64"))
    for leaf in [2, 5]:
        result.extend(
            [
                (
                    f"extra_trees_leaf{leaf}",
                    ExtraTreesClassifier(
                        n_estimators=500,
                        min_samples_leaf=leaf,
                        max_features="sqrt",
                        n_jobs=-1,
                        random_state=seed,
                    ),
                    "raw",
                ),
                (
                    f"random_forest_leaf{leaf}",
                    RandomForestClassifier(
                        n_estimators=500,
                        min_samples_leaf=leaf,
                        max_features="sqrt",
                        n_jobs=-1,
                        random_state=seed,
                    ),
                    "raw",
                ),
            ]
        )
    result.append(
        (
            "lgbm",
            lgb.LGBMClassifier(
                n_estimators=500,
                learning_rate=0.03,
                num_leaves=31,
                min_child_samples=20,
                subsample=0.8,
                colsample_bytree=0.6,
                reg_lambda=1.0,
                verbosity=-1,
                n_jobs=-1,
                random_state=seed,
            ),
            "raw",
        )
    )
    return result


def score(model, features: np.ndarray) -> np.ndarray:
    if hasattr(model, "predict_proba"):
        return model.predict_proba(features)[:, 1]
    return model.decision_function(features)


compounds = read("compounds").sort("compound_id")
all_ids = compounds.select("compound_id")
features_table = read("features")
columns = [column for column in features_table.columns if column.startswith("pxr_cp_emb_")]
features = (
    all_ids.join(features_table.select("compound_id", *columns), on="compound_id", validate="1:1")
    .select(columns)
    .to_numpy()
    .astype(np.float32)
)
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
id_to_index = {value: row for row, value in enumerate(all_ids["compound_id"])}
test_index = np.asarray([id_to_index[value] for value in test["compound_id"]])
labels = read("labels").filter(pl.col("endpoint") == cfg["endpoint"]).select("compound_id", "value")

metric_rows = []
oof_rows = []
test_rows = []
for strategy in cfg["strategies"]:
    assignment = (
        read("assignments")
        .filter(pl.col("strategy") == strategy)
        .select("compound_id", "fold", "group")
    )
    data = labels.join(assignment, on="compound_id", validate="1:1").sort("compound_id")
    data_index = np.asarray([id_to_index[value] for value in data["compound_id"]])
    y = data["value"].to_numpy().astype(int)
    folds = data["fold"].to_numpy().astype(int)
    names = [name for name, _, _ in model_specs(cfg["random_seed"])]
    predictions = {name: np.empty(len(data)) for name in names}
    test_predictions = {name: [] for name in names}
    for fold in np.unique(folds):
        train_rows = np.flatnonzero(folds != fold)
        valid_rows = np.flatnonzero(folds == fold)
        x_train = features[data_index[train_rows]]
        x_valid = features[data_index[valid_rows]]
        x_test = features[test_index]
        scaler = StandardScaler().fit(x_train)
        scaled_train = scaler.transform(x_train)
        scaled_valid = scaler.transform(x_valid)
        scaled_test = scaler.transform(x_test)
        pca = PCA(n_components=64, whiten=True, random_state=cfg["random_seed"] + fold).fit(
            scaled_train
        )
        views = {
            "raw": (x_train, x_valid, x_test),
            "scaled": (scaled_train, scaled_valid, scaled_test),
            "pca64": (
                pca.transform(scaled_train),
                pca.transform(scaled_valid),
                pca.transform(scaled_test),
            ),
        }
        for name, model, view in model_specs(cfg["random_seed"] + fold):
            train_features, valid_features, test_features = views[view]
            model.fit(train_features, y[train_rows])
            predictions[name][valid_rows] = score(model, valid_features)
            test_predictions[name].append(score(model, test_features))
        print(f"{strategy} fold {fold} complete", flush=True)

    for name, values in predictions.items():
        selected = fixed_calls(values, folds)
        metric_rows.append(
            {
                "strategy": strategy,
                "method": name,
                "AUROC": float(roc_auc_score(y, values)),
                "AP": float(average_precision_score(y, values)),
                "MCC_observed": float(matthews_corrcoef(y, selected)),
                "MCC_target_prior": expected_mcc(y, selected),
            }
        )
        averaged_test = np.mean(test_predictions[name], axis=0)
        for row, prediction in enumerate(values):
            oof_rows.append(
                {
                    "compound_id": data["compound_id"][row],
                    "endpoint": cfg["endpoint"],
                    "strategy": strategy,
                    "fold": int(folds[row]),
                    "method": name,
                    "prediction": float(prediction),
                    "y_true": int(y[row]),
                }
            )
        for row, prediction in enumerate(averaged_test):
            test_rows.append(
                {
                    "compound_id": test["compound_id"][row],
                    "endpoint": cfg["endpoint"],
                    "strategy": strategy,
                    "fold": -1,
                    "method": name,
                    "prediction": float(prediction),
                }
            )

metrics = pl.DataFrame(metric_rows).sort("strategy", "MCC_target_prior", descending=[False, True])
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(oof_rows).write_parquet(out / "oof.parquet")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "embedding_dimensions": len(columns),
            "methods": [name for name, _, _ in model_specs(cfg["random_seed"])],
            "selection": "Diagnostic model-family comparison; no submission selected.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
