"""Evaluate frozen parent/product CheMeleon embedding differences for TDI."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.decomposition import PCA
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
embeddings = read("embeddings").sort("compound_id")
scalars = read("scalars").sort("compound_id")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
if not embeddings.select("compound_id").equals(compounds.select("compound_id")):
    raise ValueError("Embedding identities differ from compounds")
if not scalars.select("compound_id").equals(compounds.select("compound_id")):
    raise ValueError("Scalar identities differ from compounds")
id_to_index = {compound_id: index for index, compound_id in enumerate(compounds["compound_id"])}
test_indices = np.asarray([id_to_index[compound_id] for compound_id in test["compound_id"]])
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in compounds["standardized_smiles"]
    ]
).astype(np.float32)


def embedding_block(scope: str) -> np.ndarray:
    columns = sorted(
        column for column in embeddings.columns if column.startswith(f"chemeleon_{scope}_")
    )
    if len(columns) != 2048:
        raise ValueError(f"Expected 2048 {scope} embedding columns")
    return embeddings.select(columns).to_numpy().astype(np.float32)


parent = embedding_block("parent")
all_products = embedding_block("all")
scope_embeddings = {"d6": embedding_block("d6"), "a4": embedding_block("a4")}
scalar_matrix = scalars.drop("compound_id").to_numpy().astype(np.float32)

metric_rows = []
oof_rows = []
test_rows = []
for endpoint in cfg["endpoints"]:
    scope = {"CYP2D6": "d6", "CYP3A4": "a4"}[endpoint.split("_", 1)[0]]
    endpoint_delta = scope_embeddings[scope] - parent
    all_delta = all_products - parent
    raw_views = {
        "parent": parent,
        "delta_endpoint": endpoint_delta,
        "delta_all_endpoint": np.column_stack([all_delta, endpoint_delta]),
        "parent_delta_endpoint": np.column_stack([parent, endpoint_delta]),
    }
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
        indices = np.asarray([id_to_index[compound_id] for compound_id in data["compound_id"]])
        folds = data["fold"].to_numpy().astype(int)
        y = data["value"].to_numpy().astype(int)
        predictions = {method: np.empty(data.height) for method in cfg["methods"]}
        test_predictions = {method: [] for method in cfg["methods"]}
        for fold in range(5):
            train_rows = np.flatnonzero(folds != fold)
            heldout_rows = np.flatnonzero(folds == fold)
            transformed = {}
            for name, matrix in raw_views.items():
                pca = PCA(
                    n_components=cfg["pca_components"],
                    svd_solver="randomized",
                    iterated_power=3,
                    random_state=cfg["seed"] + fold,
                )
                pca.fit(matrix[indices[train_rows]])
                transformed[name] = (
                    pca.transform(matrix[indices[train_rows]]).astype(np.float32),
                    pca.transform(matrix[indices[heldout_rows]]).astype(np.float32),
                    pca.transform(matrix[test_indices]).astype(np.float32),
                )
            matrices = {
                "morgan_lgbm": (
                    morgan[indices[train_rows]],
                    morgan[indices[heldout_rows]],
                    morgan[test_indices],
                ),
                "morgan_parent_pca_lgbm": tuple(
                    np.column_stack([base, extra])
                    for base, extra in zip(
                        (
                            morgan[indices[train_rows]],
                            morgan[indices[heldout_rows]],
                            morgan[test_indices],
                        ),
                        transformed["parent"],
                        strict=True,
                    )
                ),
                "morgan_product_delta_endpoint_pca_lgbm": tuple(
                    np.column_stack([base, extra])
                    for base, extra in zip(
                        (
                            morgan[indices[train_rows]],
                            morgan[indices[heldout_rows]],
                            morgan[test_indices],
                        ),
                        transformed["delta_endpoint"],
                        strict=True,
                    )
                ),
                "morgan_product_delta_all_endpoint_pca_lgbm": tuple(
                    np.column_stack([base, extra])
                    for base, extra in zip(
                        (
                            morgan[indices[train_rows]],
                            morgan[indices[heldout_rows]],
                            morgan[test_indices],
                        ),
                        transformed["delta_all_endpoint"],
                        strict=True,
                    )
                ),
                "morgan_parent_delta_endpoint_pca_lgbm": tuple(
                    np.column_stack([base, extra])
                    for base, extra in zip(
                        (
                            morgan[indices[train_rows]],
                            morgan[indices[heldout_rows]],
                            morgan[test_indices],
                        ),
                        transformed["parent_delta_endpoint"],
                        strict=True,
                    )
                ),
                "morgan_product_embedding_scalars_lgbm": (
                    np.column_stack(
                        [morgan[indices[train_rows]], scalar_matrix[indices[train_rows]]]
                    ),
                    np.column_stack(
                        [morgan[indices[heldout_rows]], scalar_matrix[indices[heldout_rows]]]
                    ),
                    np.column_stack([morgan[test_indices], scalar_matrix[test_indices]]),
                ),
            }
            for method in cfg["methods"]:
                x_train, x_heldout, x_test = matrices[method]
                model = lgb.LGBMClassifier(
                    **cfg["lgbm"], verbosity=-1, n_jobs=1, random_state=cfg["seed"] + fold
                )
                model.fit(x_train, y[train_rows])
                predictions[method][heldout_rows] = model.predict_proba(x_heldout)[:, 1]
                test_predictions[method].append(model.predict_proba(x_test)[:, 1])
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
            "pca_components": cfg["pca_components"],
            "pca_outer_training_only": True,
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
