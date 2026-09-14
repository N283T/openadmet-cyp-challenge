"""Test unhashed Morgan environments centered on TDI-alert neighborhoods."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import rankdata
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

from cyp_challenge.tdi_alerts import compiled_alerts

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = unit_rank(values[mask])
    return result


def fixed_rate_calls(scores: np.ndarray, folds: np.ndarray, rate: float) -> np.ndarray:
    calls = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = max(1, int(round(rate * len(positions))))
        calls[positions[np.argsort(-scores[positions], kind="stable")[:count]]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray, prevalence: float) -> float:
    positive = y == 1
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def fit_lr(train_dicts: list[dict], y: np.ndarray, query_dicts: list[dict]) -> np.ndarray:
    vectorizer = DictVectorizer(sparse=True)
    x_train = vectorizer.fit_transform(train_dicts)
    x_query = vectorizer.transform(query_dicts)
    for matrix in (x_train, x_query):
        matrix.indices = matrix.indices.astype(np.int32, copy=False)
        matrix.indptr = matrix.indptr.astype(np.int32, copy=False)
    model = LogisticRegression(
        C=cfg["lr_c"],
        solver="liblinear",
        class_weight="balanced",
        max_iter=4000,
        random_state=20260913,
    ).fit(x_train, y)
    return model.predict_proba(x_query)[:, 1]


compounds = read("compounds").sort("compound_id")
all_ids = compounds["compound_id"].to_list()
id_to_index = {value: row for row, value in enumerate(all_ids)}
molecules = [Chem.MolFromSmiles(value) for value in compounds["standardized_smiles"]]
alerts = compiled_alerts()
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2)
environment_rows = []
coverage = []
for compound_id, molecule in zip(all_ids, molecules, strict=True):
    alert_atoms = {
        atom
        for _, query in alerts
        for match in molecule.GetSubstructMatches(query)
        for atom in match
    }
    neighborhood_atoms = set(alert_atoms)
    for atom in alert_atoms:
        neighborhood_atoms.update(
            neighbor.GetIdx() for neighbor in molecule.GetAtomWithIdx(atom).GetNeighbors()
        )
    additional = rdFingerprintGenerator.AdditionalOutput()
    additional.AllocateBitInfoMap()
    generator.GetSparseCountFingerprint(molecule, additionalOutput=additional)
    views = {"alert_center": {}, "alert_neighborhood": {}}
    for environment_id, occurrences in additional.GetBitInfoMap().items():
        for center, radius in occurrences:
            if center in alert_atoms:
                key = f"r{radius}:{environment_id}"
                views["alert_center"][key] = views["alert_center"].get(key, 0.0) + 1.0
            if center in neighborhood_atoms:
                key = f"r{radius}:{environment_id}"
                views["alert_neighborhood"][key] = views["alert_neighborhood"].get(key, 0.0) + 1.0
    environment_rows.append(views)
    coverage.append(
        {
            "compound_id": compound_id,
            "alert_atoms": len(alert_atoms),
            "neighborhood_atoms": len(neighborhood_atoms),
        }
    )
pl.DataFrame(coverage).write_parquet(out / "coverage.parquet")

block_oof = read("block_oof")
block_test = read("block_test")
metric_rows, oof_rows, test_rows = [], [], []
for endpoint in cfg["endpoints"]:
    settings = cfg["endpoint_settings"][endpoint]
    allowed_radii = set(settings["radii"])
    for strategy in cfg["strategies"]:
        base = block_oof.filter(
            (pl.col("endpoint") == endpoint)
            & (pl.col("strategy") == strategy)
            & (pl.col("method") == settings["block_method"])
        ).sort("compound_id")
        deploy = block_test.filter(
            (pl.col("endpoint") == endpoint)
            & (pl.col("strategy") == strategy)
            & (pl.col("method") == settings["block_method"])
        ).sort("compound_id")
        folds = base["fold"].to_numpy().astype(int)
        y = base["y_true"].to_numpy().astype(int)
        train_indices = np.asarray([id_to_index[value] for value in base["compound_id"]])
        test_indices = np.asarray([id_to_index[value] for value in deploy["compound_id"]])
        base_rank = fold_rank(base["prediction"].to_numpy(), folds)
        base_test_rank = unit_rank(deploy["prediction"].to_numpy())
        candidates = {"challenge_block": (base_rank, base_test_rank, 1, 0)}
        for view in ["alert_center", "alert_neighborhood"]:
            dictionaries = [
                {
                    key: value
                    for key, value in row[view].items()
                    if int(key.split(":", 1)[0][1:]) in allowed_radii
                }
                for row in environment_rows
            ]
            oof_prediction = np.empty(len(base))
            fold_test = []
            for fold in range(5):
                train = np.flatnonzero(folds != fold)
                heldout = np.flatnonzero(folds == fold)
                oof_prediction[heldout] = fit_lr(
                    [dictionaries[train_indices[row]] for row in train],
                    y[train],
                    [dictionaries[train_indices[row]] for row in heldout],
                )
                fold_test.append(
                    fit_lr(
                        [dictionaries[train_indices[row]] for row in train],
                        y[train],
                        [dictionaries[row] for row in test_indices],
                    )
                )
            local_rank = fold_rank(oof_prediction, folds)
            local_test_rank = unit_rank(np.mean(fold_test, axis=0))
            candidates[view] = (local_rank, local_test_rank, 0, 1)
            for base_weight, local_weight in cfg["blend_weights"]:
                candidates[f"challenge{base_weight}_{view}{local_weight}"] = (
                    (base_weight * base_rank + local_weight * local_rank)
                    / (base_weight + local_weight),
                    (base_weight * base_test_rank + local_weight * local_test_rank)
                    / (base_weight + local_weight),
                    base_weight,
                    local_weight,
                )
        for method, (
            prediction,
            deploy_prediction,
            base_weight,
            local_weight,
        ) in candidates.items():
            calls = fixed_rate_calls(prediction, folds, settings["call_rate"])
            metric_rows.append(
                {
                    "endpoint": endpoint,
                    "strategy": strategy,
                    "method": method,
                    "base_weight": base_weight,
                    "local_weight": local_weight,
                    "AUROC": float(roc_auc_score(y, prediction)),
                    "AP": float(average_precision_score(y, prediction)),
                    "MCC_observed": float(matthews_corrcoef(y, calls)),
                    "MCC_target_prior": expected_mcc(y, calls, settings["target_prior"]),
                }
            )
            for row, value in enumerate(prediction):
                oof_rows.append(
                    {
                        "compound_id": base["compound_id"][row],
                        "endpoint": endpoint,
                        "strategy": strategy,
                        "fold": int(folds[row]),
                        "group": base["group"][row],
                        "y_true": int(y[row]),
                        "method": method,
                        "prediction": float(value),
                        "call": bool(calls[row]),
                    }
                )
            for row, value in enumerate(deploy_prediction):
                test_rows.append(
                    {
                        "compound_id": deploy["compound_id"][row],
                        "endpoint": endpoint,
                        "strategy": strategy,
                        "method": method,
                        "prediction": float(value),
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
            "scope": "Alert-centered exact Morgan environments; dual-split exploratory screen.",
            "selection": "No candidate packaging or submission.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
