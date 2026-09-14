"""Evaluate an independent FDA CYP3A4 TDI task as an A4 prediction block."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import rankdata
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from sklearn.svm import SVC

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


def dense_bits(smiles: list[str], generator) -> tuple[np.ndarray, list]:
    matrix = np.zeros((len(smiles), cfg["fp_size"]), dtype=np.float32)
    fingerprints = []
    for row, value in enumerate(smiles):
        molecule = Chem.MolFromSmiles(value)
        fingerprint = generator.GetFingerprint(molecule)
        DataStructs.ConvertToNumpyArray(fingerprint, matrix[row])
        fingerprints.append(fingerprint)
    return matrix, fingerprints


def tanimoto_matrix(query: list, reference: list) -> np.ndarray:
    return np.asarray(
        [DataStructs.BulkTanimotoSimilarity(fp, reference) for fp in query], dtype=np.float32
    )


external = read("external_records").filter(
    pl.col("value").is_in([0, 1]) & pl.col("challenge_compound_id").is_null()
)
compounds = read("compounds").sort("compound_id")
generator = rdFingerprintGenerator.GetMorganGenerator(radius=cfg["radius"], fpSize=cfg["fp_size"])
x_external, fp_external = dense_bits(external["canonical_smiles"].to_list(), generator)
x_challenge, fp_challenge = dense_bits(compounds["standardized_smiles"].to_list(), generator)
y_external = external["value"].to_numpy().astype(int)

models = {
    "fda_lr": LogisticRegression(
        C=cfg["lr_c"],
        solver="liblinear",
        class_weight="balanced",
        max_iter=4000,
        random_state=20260913,
    ).fit(x_external, y_external),
    "fda_extra_trees": ExtraTreesClassifier(
        n_estimators=cfg["trees"],
        min_samples_leaf=cfg["min_samples_leaf"],
        max_features="sqrt",
        class_weight="balanced",
        n_jobs=-1,
        random_state=20260913,
    ).fit(x_external, y_external),
}
external_predictions = {
    "fda_lr": models["fda_lr"].predict_proba(x_challenge)[:, 1],
    "fda_extra_trees": models["fda_extra_trees"].predict_proba(x_challenge)[:, 1],
}
kernel_train = tanimoto_matrix(fp_external, fp_external)
kernel_query = tanimoto_matrix(fp_challenge, fp_external)
svm = SVC(C=cfg["svm_c"], kernel="precomputed", class_weight="balanced").fit(
    kernel_train, y_external
)
external_predictions["fda_tanimoto_svm"] = svm.decision_function(kernel_query)
id_to_index = {value: row for row, value in enumerate(compounds["compound_id"])}

block_oof = read("block_oof").filter(
    (pl.col("endpoint") == cfg["endpoint"]) & (pl.col("method") == cfg["block_method"])
)
block_test = read("block_test").filter(
    (pl.col("endpoint") == cfg["endpoint"]) & (pl.col("method") == cfg["block_method"])
)
metric_rows = []
oof_rows = []
test_rows = []
for strategy in cfg["strategies"]:
    base = block_oof.filter(pl.col("strategy") == strategy).sort("compound_id")
    deploy = block_test.filter(pl.col("strategy") == strategy).sort("compound_id")
    folds = base["fold"].to_numpy().astype(int)
    y = base["y_true"].to_numpy().astype(int)
    base_rank = fold_rank(base["prediction"].to_numpy(), folds)
    base_test_rank = unit_rank(deploy["prediction"].to_numpy())
    oof_indices = np.asarray([id_to_index[value] for value in base["compound_id"]])
    test_indices = np.asarray([id_to_index[value] for value in deploy["compound_id"]])
    candidates = {"challenge_block": (base_rank, base_test_rank, 1, 0)}
    for external_name, prediction in external_predictions.items():
        external_rank = fold_rank(prediction[oof_indices], folds)
        external_test_rank = unit_rank(prediction[test_indices])
        candidates[external_name] = (external_rank, external_test_rank, 0, 1)
        for base_weight, external_weight in cfg["blend_weights"]:
            name = f"challenge{base_weight}_{external_name}{external_weight}"
            candidates[name] = (
                (base_weight * base_rank + external_weight * external_rank)
                / (base_weight + external_weight),
                (base_weight * base_test_rank + external_weight * external_test_rank)
                / (base_weight + external_weight),
                base_weight,
                external_weight,
            )
    for method, (prediction, deploy_prediction, base_weight, external_weight) in candidates.items():
        calls = fixed_rate_calls(prediction, folds, cfg["call_rate"])
        metric_rows.append(
            {
                "strategy": strategy,
                "method": method,
                "base_weight": base_weight,
                "external_weight": external_weight,
                "AUROC": float(roc_auc_score(y, prediction)),
                "AP": float(average_precision_score(y, prediction)),
                "MCC_observed": float(matthews_corrcoef(y, calls)),
                "MCC_target_prior": expected_mcc(y, calls, cfg["target_prior"]),
            }
        )
        for row, value in enumerate(prediction):
            oof_rows.append(
                {
                    "compound_id": base["compound_id"][row],
                    "strategy": strategy,
                    "endpoint": cfg["endpoint"],
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
                    "strategy": strategy,
                    "method": method,
                    "prediction": float(value),
                }
            )

metrics = pl.DataFrame(metric_rows).sort("strategy", "MCC_target_prior", descending=[False, True])
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(oof_rows).write_parquet(out / "oof.parquet")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "external_training_rows": external.height,
            "excluded_exact_challenge_overlaps": 40,
            "scope": (
                "Independent-task prediction block; no source labels merged into "
                "challenge training."
            ),
            "selection": "Exploratory only; no candidate packaging or submission.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
