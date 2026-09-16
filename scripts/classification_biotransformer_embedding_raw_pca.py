"""Compare raw and fold-local PCA CheMeleon parent/product embeddings."""

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
    result = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(rate * len(rows))))
        result[rows[np.argsort(-scores[rows], kind="stable")[:count]]] = True
    return result


def expected_mcc(y: np.ndarray, selected: np.ndarray, prevalence: float) -> float:
    positive = y == 1
    tpr, fpr = float(selected[positive].mean()), float(selected[~positive].mean())
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


compounds = read("compounds").sort("compound_id")
labels = read("labels")
assignments = read("assignments")
embeddings = read("embeddings").sort("compound_id")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
if not embeddings.select("compound_id").equals(compounds.select("compound_id")):
    raise ValueError("Embedding identities differ from compounds")
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
        raise ValueError(f"Expected 2048 {scope} columns, got {len(columns)}")
    return embeddings.select(columns).to_numpy().astype(np.float32)


parent = embedding_block("parent")
scope_embeddings = {"d6": embedding_block("d6"), "a4": embedding_block("a4")}
methods = ["morgan_lgbm"]
for view in cfg["views"]:
    methods.append(f"morgan_{view}_raw_lgbm")
    methods.extend(
        f"morgan_{view}_pca{components}_lgbm"
        for components in cfg["pca_components"]
    )

metrics = []
oof_rows = []
test_rows = []
for endpoint in cfg["endpoints"]:
    scope = {"CYP2D6": "d6", "CYP3A4": "a4"}[endpoint.split("_", 1)[0]]
    delta = scope_embeddings[scope] - parent
    raw_views = {
        "parent": parent,
        "delta_endpoint": delta,
        "parent_delta_endpoint": np.column_stack([parent, delta]),
    }
    endpoint_labels = labels.filter(pl.col("endpoint") == endpoint).select(
        "compound_id", "value"
    )
    settings = cfg["endpoint_settings"][endpoint]
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
        predictions = {method: np.empty(data.height) for method in methods}
        test_predictions = {method: [] for method in methods}
        for fold in range(5):
            train_rows = np.flatnonzero(folds != fold)
            heldout_rows = np.flatnonzero(folds == fold)
            matrices = {
                "morgan_lgbm": (
                    morgan[indices[train_rows]],
                    morgan[indices[heldout_rows]],
                    morgan[test_indices],
                )
            }
            for view, raw in raw_views.items():
                matrices[f"morgan_{view}_raw_lgbm"] = tuple(
                    np.column_stack([base, extra])
                    for base, extra in zip(
                        (
                            morgan[indices[train_rows]],
                            morgan[indices[heldout_rows]],
                            morgan[test_indices],
                        ),
                        (
                            raw[indices[train_rows]],
                            raw[indices[heldout_rows]],
                            raw[test_indices],
                        ),
                        strict=True,
                    )
                )
                for components in cfg["pca_components"]:
                    pca = PCA(
                        n_components=components,
                        svd_solver="randomized",
                        iterated_power=3,
                        random_state=cfg["seed"] + fold,
                    )
                    pca.fit(raw[indices[train_rows]])
                    matrices[f"morgan_{view}_pca{components}_lgbm"] = tuple(
                        np.column_stack([base, pca.transform(extra).astype(np.float32)])
                        for base, extra in zip(
                            (
                                morgan[indices[train_rows]],
                                morgan[indices[heldout_rows]],
                                morgan[test_indices],
                            ),
                            (
                                raw[indices[train_rows]],
                                raw[indices[heldout_rows]],
                                raw[test_indices],
                            ),
                            strict=True,
                        )
                    )
            for method in methods:
                x_train, x_heldout, x_test = matrices[method]
                model = lgb.LGBMClassifier(
                    **cfg["lgbm"],
                    verbosity=-1,
                    n_jobs=1,
                    random_state=cfg["seed"] + fold,
                )
                model.fit(x_train, y[train_rows])
                predictions[method][heldout_rows] = model.predict_proba(x_heldout)[:, 1]
                test_predictions[method].append(model.predict_proba(x_test)[:, 1])
            print(f"{endpoint} {strategy} fold {fold} complete", flush=True)
        for method in methods:
            probability = predictions[method]
            calls = fixed_rate_calls(probability, folds, settings["call_rate"])
            metrics.append(
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
            mean_test = np.mean(test_predictions[method], axis=0)
            for row, score in enumerate(mean_test):
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

metric_table = pl.DataFrame(metrics).sort("endpoint", "strategy", "method")
metric_table.write_csv(out / "metrics.csv")
pl.DataFrame(oof_rows).write_parquet(out / "oof.parquet")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "selection": cfg["selection"],
            "methods": methods,
            "pca_outer_training_only": True,
            "raw_dimensions": {
                "parent": 2048,
                "delta_endpoint": 2048,
                "parent_delta_endpoint": 4096,
            },
        },
        indent=2,
    )
    + "\n"
)
print(metric_table, flush=True)
