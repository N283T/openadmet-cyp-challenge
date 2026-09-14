"""Compact classification screen on UMAP and scaffold folds."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.special import expit, logit
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from sklearn.preprocessing import StandardScaler

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def aligned_matrix(ids: pl.DataFrame, table: pl.DataFrame) -> np.ndarray:
    columns = [column for column in table.columns if column.startswith("emb_")]
    values = (
        ids.join(
            table.select("compound_id", *columns),
            on="compound_id",
            how="left",
            validate="1:1",
            maintain_order="left",
        )
        .select(columns)
        .to_numpy()
        .astype(np.float32)
    )
    if not np.isfinite(values).all():
        raise ValueError("Embedding table does not fully cover compound IDs")
    return values


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


def fit_predict(
    method: str,
    features: dict[str, np.ndarray],
    train: np.ndarray,
    target: np.ndarray,
    query: np.ndarray,
) -> np.ndarray:
    if method == "morgan_r2_lr":
        model = LogisticRegression(
            C=0.05, solver="liblinear", max_iter=4000, random_state=20260913
        ).fit(features["morgan"][train], target)
        return model.predict_proba(features["morgan"][query])[:, 1]

    view = "class" if "class" in method else "mixed"
    scaler = StandardScaler().fit(features[view][train])
    x_train = scaler.transform(features[view][train])
    x_query = scaler.transform(features[view][query])
    if method.startswith("pls16"):
        pls = PLSRegression(n_components=16, scale=False, max_iter=1000).fit(x_train, target)
        x_train = pls.transform(x_train)
        x_query = pls.transform(x_query)
        c_value = 1.0
    else:
        c_value = 0.05
    model = LogisticRegression(
        C=c_value, solver="liblinear", max_iter=4000, random_state=20260913
    ).fit(x_train, target)
    return model.predict_proba(x_query)[:, 1]


def residual_adjustment(
    query_indices: np.ndarray,
    pool_indices: np.ndarray,
    query_probability: np.ndarray,
    residuals: np.ndarray,
    fingerprints: list,
) -> np.ndarray:
    pool_fingerprints = [fingerprints[index] for index in pool_indices]
    corrections = np.empty(len(query_indices), dtype=float)
    for row, index in enumerate(query_indices):
        similarities = np.asarray(
            DataStructs.BulkTanimotoSimilarity(fingerprints[index], pool_fingerprints)
        )
        nearest = np.argsort(-similarities, kind="stable")[: cfg["residual_k"]]
        weights = similarities[nearest] ** cfg["residual_similarity_power"]
        correction = (
            float(np.average(residuals[nearest], weights=weights)) if weights.sum() > 0 else 0.0
        )
        corrections[row] = similarities[nearest[0]] * correction
    clipped = np.clip(query_probability, 1e-6, 1 - 1e-6)
    return expit(logit(clipped) + cfg["residual_lambda"] * corrections)


compounds = read("compounds").sort("compound_id")
all_ids = compounds.select("compound_id")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
labels = read("labels")
assignments = read("assignments")
id_to_index = {compound_id: row for row, compound_id in enumerate(all_ids["compound_id"])}
test_index = np.asarray([id_to_index[value] for value in test["compound_id"]])
molecules = [Chem.MolFromSmiles(value) for value in compounds["standardized_smiles"]]
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fingerprints = [generator.GetFingerprint(molecule) for molecule in molecules]
features = {
    "morgan": np.stack(
        [generator.GetFingerprintAsNumPy(molecule) for molecule in molecules]
    ).astype(np.float32),
    "class": aligned_matrix(all_ids, read("class_emb")),
    "mixed": aligned_matrix(all_ids, read("mixed_emb")),
}

oof_rows = []
test_rows = []
metric_rows = []
for endpoint in cfg["endpoints"]:
    endpoint_cfg = cfg["endpoint_settings"][endpoint]
    endpoint_labels = labels.filter(pl.col("endpoint") == endpoint).select("compound_id", "value")
    for strategy in cfg["strategies"]:
        fold_table = assignments.filter(pl.col("strategy") == strategy).select(
            "compound_id", "fold", "group"
        )
        data = endpoint_labels.join(fold_table, on="compound_id", validate="1:1").sort(
            "compound_id"
        )
        data_index = np.asarray([id_to_index[value] for value in data["compound_id"]])
        y = data["value"].to_numpy().astype(int)
        folds = data["fold"].to_numpy().astype(int)
        predictions = {method: np.empty(len(data), dtype=float) for method in cfg["methods"]}
        test_predictions = {method: [] for method in cfg["methods"]}

        for outer in range(5):
            train_rows = np.flatnonzero(folds != outer)
            heldout_rows = np.flatnonzero(folds == outer)
            for method in cfg["methods"]:
                val_probability = fit_predict(
                    method,
                    features,
                    data_index[train_rows],
                    y[train_rows],
                    data_index[heldout_rows],
                )
                fold_test_probability = fit_predict(
                    method,
                    features,
                    data_index[train_rows],
                    y[train_rows],
                    test_index,
                )
                predictions[method][heldout_rows] = val_probability
                test_predictions[method].append(fold_test_probability)

                if method not in cfg["residual_base_methods"]:
                    continue
                inner_probability = np.empty(len(train_rows), dtype=float)
                for inner in sorted(set(folds[train_rows])):
                    inner_heldout = np.flatnonzero(folds[train_rows] == inner)
                    inner_train = np.flatnonzero(folds[train_rows] != inner)
                    inner_probability[inner_heldout] = fit_predict(
                        method,
                        features,
                        data_index[train_rows[inner_train]],
                        y[train_rows[inner_train]],
                        data_index[train_rows[inner_heldout]],
                    )
                residuals = y[train_rows] - inner_probability
                adjusted_name = f"{method}_local_residual"
                predictions.setdefault(adjusted_name, np.empty(len(data)))[heldout_rows] = (
                    residual_adjustment(
                        data_index[heldout_rows],
                        data_index[train_rows],
                        val_probability,
                        residuals,
                        fingerprints,
                    )
                )
                test_predictions.setdefault(adjusted_name, []).append(
                    residual_adjustment(
                        test_index,
                        data_index[train_rows],
                        fold_test_probability,
                        residuals,
                        fingerprints,
                    )
                )
            print(f"{endpoint} {strategy} fold {outer} complete", flush=True)

        for method, probability in predictions.items():
            calls = fixed_rate_calls(probability, folds, endpoint_cfg["call_rate"])
            metric_rows.append(
                {
                    "endpoint": endpoint,
                    "strategy": strategy,
                    "method": method,
                    "AUROC": float(roc_auc_score(y, probability)),
                    "AP": float(average_precision_score(y, probability)),
                    "MCC_observed": float(matthews_corrcoef(y, calls)),
                    "MCC_target_prior": expected_mcc(y, calls, endpoint_cfg["target_prior"]),
                    "calls": int(calls.sum()),
                }
            )
            for row, prediction in enumerate(probability):
                oof_rows.append(
                    {
                        "compound_id": data["compound_id"][row],
                        "endpoint": endpoint,
                        "strategy": strategy,
                        "fold": int(folds[row]),
                        "group": data["group"][row],
                        "y_true": int(y[row]),
                        "method": method,
                        "prediction": float(prediction),
                        "call": bool(calls[row]),
                    }
                )
            averaged_test = np.mean(test_predictions[method], axis=0)
            for row, prediction in enumerate(averaged_test):
                test_rows.append(
                    {
                        "compound_id": test["compound_id"][row],
                        "endpoint": endpoint,
                        "strategy": strategy,
                        "method": method,
                        "prediction": float(prediction),
                        "split_id": f"{strategy}_five_fold",
                        "fold": -1,
                    }
                )

metrics = pl.DataFrame(metric_rows).sort(
    "endpoint", "strategy", "MCC_target_prior", descending=[False, False, True]
)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(oof_rows).write_parquet(out / "oof.parquet")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "strategies": cfg["strategies"],
            "methods": cfg["methods"],
            "local_residual": {
                "base_methods": cfg["residual_base_methods"],
                "k": cfg["residual_k"],
                "similarity_power": cfg["residual_similarity_power"],
                "lambda": cfg["residual_lambda"],
                "honesty": (
                    "Outer-training residuals come from inner models excluding both "
                    "the outer fold and each residual row's inner fold."
                ),
            },
            "selection": "Prespecified comparison; no method is selected for submission.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics.write_csv(), flush=True)
