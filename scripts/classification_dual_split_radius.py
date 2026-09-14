"""Endpoint-specific Morgan radius and unhashed-environment comparison."""

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


def fit_lr(x_train, y_train, x_query):
    model = LogisticRegression(
        C=0.05, solver="liblinear", max_iter=4000, random_state=20260913
    ).fit(x_train, y_train)
    return model.predict_proba(x_query)[:, 1]


compounds = read("compounds").sort("compound_id")
all_ids = compounds.select("compound_id")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
labels = read("labels")
assignments = read("assignments")
id_to_index = {compound_id: row for row, compound_id in enumerate(all_ids["compound_id"])}
test_index = np.asarray([id_to_index[value] for value in test["compound_id"]])
molecules = [Chem.MolFromSmiles(value) for value in compounds["standardized_smiles"]]
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
all_bits = np.zeros((len(molecules), 2048), dtype=np.float32)
radius_bits = [np.zeros_like(all_bits) for _ in range(3)]
unhashed_by_radius = [[] for _ in range(3)]
for row, molecule in enumerate(molecules):
    additional = rdFingerprintGenerator.AdditionalOutput()
    additional.AllocateBitInfoMap()
    fingerprint = generator.GetFingerprint(molecule, additionalOutput=additional)
    all_bits[row, list(fingerprint.GetOnBits())] = 1
    for bit, environments in additional.GetBitInfoMap().items():
        for _, radius in environments:
            radius_bits[radius][row, bit] = 1
    radius_dicts = [{} for _ in range(3)]
    sparse_additional = rdFingerprintGenerator.AdditionalOutput()
    sparse_additional.AllocateBitInfoMap()
    sparse = generator.GetSparseCountFingerprint(molecule, additionalOutput=sparse_additional)
    sparse_counts = sparse.GetNonzeroElements()
    for environment_id, environments in sparse_additional.GetBitInfoMap().items():
        radii = {radius for _, radius in environments}
        for radius in radii:
            radius_dicts[radius][str(environment_id)] = float(sparse_counts[environment_id])
    for radius in range(3):
        unhashed_by_radius[radius].append(radius_dicts[radius])

oof_rows = []
test_rows = []
metric_rows = []
coverage_rows = []
for endpoint in cfg["endpoints"]:
    settings = cfg["endpoint_settings"][endpoint]
    radii = settings["radii"]
    specific_bits = np.maximum.reduce([radius_bits[radius] for radius in radii])
    specific_dicts = []
    for row in range(len(molecules)):
        combined = {}
        for radius in radii:
            combined.update(
                {
                    f"r{radius}:{key}": value
                    for key, value in unhashed_by_radius[radius][row].items()
                }
            )
        specific_dicts.append(combined)
    endpoint_labels = labels.filter(pl.col("endpoint") == endpoint).select("compound_id", "value")
    for strategy in cfg["strategies"]:
        split = assignments.filter(pl.col("strategy") == strategy).select(
            "compound_id", "fold", "group"
        )
        data = endpoint_labels.join(split, on="compound_id", validate="1:1").sort("compound_id")
        data_index = np.asarray([id_to_index[value] for value in data["compound_id"]])
        y = data["value"].to_numpy().astype(int)
        folds = data["fold"].to_numpy().astype(int)
        predictions = {
            "morgan_all_bit": np.empty(len(data)),
            "radius_specific_bit": np.empty(len(data)),
            "radius_specific_unhashed_count": np.empty(len(data)),
        }
        test_predictions = {method: [] for method in predictions}
        for fold in range(5):
            train_rows = np.flatnonzero(folds != fold)
            heldout_rows = np.flatnonzero(folds == fold)
            train_index = data_index[train_rows]
            heldout_index = data_index[heldout_rows]
            predictions["morgan_all_bit"][heldout_rows] = fit_lr(
                all_bits[train_index], y[train_rows], all_bits[heldout_index]
            )
            test_predictions["morgan_all_bit"].append(
                fit_lr(all_bits[train_index], y[train_rows], all_bits[test_index])
            )
            predictions["radius_specific_bit"][heldout_rows] = fit_lr(
                specific_bits[train_index],
                y[train_rows],
                specific_bits[heldout_index],
            )
            test_predictions["radius_specific_bit"].append(
                fit_lr(specific_bits[train_index], y[train_rows], specific_bits[test_index])
            )
            vectorizer = DictVectorizer(sparse=True)
            x_train = vectorizer.fit_transform([specific_dicts[index] for index in train_index])
            x_heldout = vectorizer.transform([specific_dicts[index] for index in heldout_index])
            x_test = vectorizer.transform([specific_dicts[index] for index in test_index])
            for matrix in [x_train, x_heldout, x_test]:
                matrix.indices = matrix.indices.astype(np.int32, copy=False)
                matrix.indptr = matrix.indptr.astype(np.int32, copy=False)
            predictions["radius_specific_unhashed_count"][heldout_rows] = fit_lr(
                x_train, y[train_rows], x_heldout
            )
            test_predictions["radius_specific_unhashed_count"].append(
                fit_lr(x_train, y[train_rows], x_test)
            )
            coverage_rows.append(
                {
                    "endpoint": endpoint,
                    "strategy": strategy,
                    "fold": fold,
                    "radii": ",".join(map(str, radii)),
                    "train_unhashed_features": len(vectorizer.vocabulary_),
                    "heldout_nonzero": int(x_heldout.nnz),
                    "test_nonzero": int(x_test.nnz),
                }
            )

        predictions["rank_mean_all_specific_bit"] = (
            fold_rank(predictions["morgan_all_bit"], folds)
            + fold_rank(predictions["radius_specific_bit"], folds)
        ) / 2
        predictions["rank_mean_all_unhashed"] = (
            fold_rank(predictions["morgan_all_bit"], folds)
            + fold_rank(predictions["radius_specific_unhashed_count"], folds)
        ) / 2
        for name in ["rank_mean_all_specific_bit", "rank_mean_all_unhashed"]:
            component = name.removeprefix("rank_mean_all_")
            component_name = (
                "radius_specific_bit"
                if component == "specific_bit"
                else "radius_specific_unhashed_count"
            )
            test_predictions[name] = [
                (
                    unit_rank(test_predictions["morgan_all_bit"][fold])
                    + unit_rank(test_predictions[component_name][fold])
                )
                / 2
                for fold in range(5)
            ]

        for method, probability in predictions.items():
            calls = fixed_rate_calls(probability, folds, settings["call_rate"])
            metric_rows.append(
                {
                    "endpoint": endpoint,
                    "strategy": strategy,
                    "radii": ",".join(map(str, radii)),
                    "method": method,
                    "AUROC": float(roc_auc_score(y, probability)),
                    "AP": float(average_precision_score(y, probability)),
                    "MCC_observed": float(matthews_corrcoef(y, calls)),
                    "MCC_target_prior": expected_mcc(y, calls, settings["target_prior"]),
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
                        "fold": -1,
                    }
                )
        print(f"{endpoint} {strategy} complete", flush=True)

metrics = pl.DataFrame(metric_rows).sort(
    "endpoint", "strategy", "MCC_target_prior", descending=[False, False, True]
)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(coverage_rows).write_csv(out / "unhashed-coverage.csv")
pl.DataFrame(oof_rows).write_parquet(out / "oof.parquet")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "D6_radii": cfg["endpoint_settings"]["CYP2D6_is_TDI"]["radii"],
            "A4_radii": cfg["endpoint_settings"]["CYP3A4_is_TDI"]["radii"],
            "unhashed_vocabulary": (
                "Fit on each outer-training fold only; unseen environments ignored."
            ),
            "selection": "Five prespecified methods; no submission candidate selected.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics.write_csv(), flush=True)
