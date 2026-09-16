"""Ablate LightGBM members inside the fixed deployed D6 hierarchy."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoint = "CYP2D6_is_TDI"


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


def fixed_calls(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(rows))))
        result[rows[np.argsort(-values[rows], kind="stable")[:count]]] = True
    return result


def top_calls(values: np.ndarray) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    result[np.argsort(-values, kind="stable")[: cfg["test_calls"]]] = True
    return result


def expected_mcc(y: np.ndarray, selected: np.ndarray) -> float:
    positive = y == 1
    tpr, fpr = float(selected[positive].mean()), float(selected[~positive].mean())
    prevalence = cfg["target_prior"]
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def model_score(model: object, features: np.ndarray) -> np.ndarray:
    if hasattr(model, "predict_proba"):
        return model.predict_proba(features)[:, 1]
    return model.decision_function(features)


reference = (
    read("reference_oof")
    .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == cfg["reference_method"]))
    .sort("compound_id")
)
reference_test = (
    read("reference_test")
    .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == cfg["reference_method"]))
    .sort("compound_id")
)
ids, test_ids = reference.select("compound_id"), reference_test.select("compound_id")
folds = reference["fold"].to_numpy().astype(int)
y = reference["y_true"].to_numpy().astype(int)


def atomic(source: str, method: str, column: str, test: bool = False) -> np.ndarray:
    table = read(source).filter((pl.col("endpoint") == endpoint) & (pl.col("method") == method))
    if test:
        table = table.group_by("compound_id").agg(pl.col(column).mean().alias("prediction"))
        aligned = test_ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    else:
        table = table.select("compound_id", pl.col(column).alias("prediction"))
        aligned = ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    return aligned["prediction"].to_numpy()


specs = [
    ("lgbm::morgan_plain", "lgbm_oof", "lgbm_test", "morgan_plain", "prediction"),
    ("lgbm::mordred_plain", "lgbm_oof", "lgbm_test", "mordred_plain", "prediction"),
    ("lgbm::fp_mordred_plain", "lgbm_oof", "lgbm_test", "fp_mordred_plain", "prediction"),
    ("lgbm::external_class_plain", "lgbm_oof", "lgbm_test", "external_class_plain", "prediction"),
    ("lgbm::external_mixed_plain", "lgbm_oof", "lgbm_test", "external_mixed_plain", "prediction"),
    (
        "lgbm::external_class_mordred_plain",
        "lgbm_oof",
        "lgbm_test",
        "external_class_mordred_plain",
        "prediction",
    ),
    ("pfn::external_class_pca256", "pfn_oof", "pfn_test", "external_class_pca256", "probability"),
    ("pfn::external_mixed_pca256", "pfn_oof", "pfn_test", "external_mixed_pca256", "probability"),
    ("pfn::chemeleon_pca256", "pfn_oof", "pfn_test", "chemeleon_pca256", "probability"),
    ("graph::tdi_only", "graph_oof", "graph_test", "tdi_only", "probability"),
    ("graph::tdi_aux", "graph_oof", "graph_test", "tdi_aux", "probability"),
]
names = [spec[0] for spec in specs]
oof_atoms = {
    name: atomic(oof_source, method, column) for name, oof_source, _, method, column in specs
}
test_atoms = {
    name: atomic(test_source, method, column, True)
    for name, _, test_source, method, column in specs
}
reconstructed = np.mean(list(oof_atoms.values()), axis=0)
if not np.allclose(reconstructed, reference["prediction"].to_numpy(), atol=1e-7):
    raise ValueError("Atomic mean does not reconstruct old11 OOF")

fixed = read("fixed_oof").sort("compound_id")
fixed_oof, fixed_test = {}, {}
compounds = read("compounds").sort("compound_id")
test_mapping = read("test_mapping").select("compound_id").unique().sort("compound_id")
compound_index = {value: index for index, value in enumerate(compounds["compound_id"])}
test_index = np.asarray([compound_index[value] for value in test_mapping["compound_id"]])
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in compounds["standardized_smiles"]
    ]
).astype(np.float32)
for short, method, model_key in [
    ("logistic", "logistic_c0.05", "logistic_models"),
    ("complement", "complement_nb_a1", "complement_models"),
    ("rbf", "rbf_svm_c1", "rbf_models"),
]:
    table = fixed.filter((pl.col("endpoint") == endpoint) & (pl.col("method") == method)).sort(
        "compound_id"
    )
    fixed_oof[short] = fold_rank(table["prediction"].to_numpy(), folds)
    fold_predictions = []
    for model_file in cfg[model_key]:
        fold_predictions.append(
            unit_rank(model_score(joblib.load(inp / model_file), morgan[test_index]))
        )
    fixed_test[short] = unit_rank(np.mean(fold_predictions, axis=0))

rows, prediction_rows, test_rows = [], [], []
scores, calls_by_variant = {}, {}
for variant, removed in cfg["ablations"].items():
    retained = [name for name in names if name not in removed]
    old_oof = fold_rank(np.mean([oof_atoms[name] for name in retained], axis=0), folds)
    old_test = unit_rank(np.mean([test_atoms[name] for name in retained], axis=0))
    live_oof = (2 * old_oof + fixed_oof["logistic"]) / 3
    live_test = (2 * old_test + fixed_test["logistic"]) / 3
    tail_oof = {
        "rbf": (2 * live_oof + fixed_oof["rbf"]) / 3,
        "complement": (2 * live_oof + fixed_oof["complement"]) / 3,
    }
    score = np.empty(len(y), dtype=float)
    for fold, tail in enumerate(cfg["current_fold_tail"]):
        mask = folds == fold
        score[mask] = tail_oof[tail][mask]
    test_score = (2 * live_test + fixed_test["complement"]) / 3
    selected, test_selected = fixed_calls(score, folds), top_calls(test_score)
    scores[variant], calls_by_variant[variant] = score, selected
    rows.append(
        {
            "variant": variant,
            "old11_members": len(retained),
            "AUROC": float(roc_auc_score(y, score)),
            "AP": float(average_precision_score(y, score)),
            "MCC_observed": float(matthews_corrcoef(y, selected)),
            "MCC_target_prior": expected_mcc(y, selected),
            "test_spearman_to_full11": float("nan") if variant == "full11" else 0.0,
            "test_call_swaps_each_way": int(
                (test_selected & ~top_calls(reference_test["prediction"].to_numpy())).sum()
            ),
        }
    )
    for index, compound_id in enumerate(ids["compound_id"]):
        prediction_rows.append(
            {
                "compound_id": compound_id,
                "fold": int(folds[index]),
                "group": fixed.filter(pl.col("method") == "logistic_c0.05").sort("compound_id")[
                    "group"
                ][index],
                "y_true": int(y[index]),
                "variant": variant,
                "prediction": float(score[index]),
                "call": bool(selected[index]),
            }
        )
    test_rows.append((variant, test_score, test_selected))

full_test = next(score for variant, score, _ in test_rows if variant == "full11")
full_test_calls = top_calls(full_test)
for row, (_variant, test_score, test_selected) in zip(rows, test_rows, strict=True):
    row["test_spearman_to_full11"] = float(spearmanr(test_score, full_test).statistic)
    row["test_call_swaps_each_way"] = int((test_selected & ~full_test_calls).sum())

groups = np.asarray([row["group"] for row in prediction_rows[: len(y)]])
group_values = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
rng = np.random.default_rng(cfg["bootstrap_seed"])
bootstrap_rows = []
for variant in cfg["ablations"]:
    if variant == "full11":
        continue
    deltas = []
    for _ in range(cfg["bootstrap_replicates"]):
        sampled = rng.choice(group_values, size=len(group_values), replace=True)
        sampled_rows = np.concatenate([group_rows[group] for group in sampled])
        deltas.append(
            expected_mcc(y[sampled_rows], calls_by_variant[variant][sampled_rows])
            - expected_mcc(y[sampled_rows], calls_by_variant["full11"][sampled_rows])
        )
    values = np.asarray(deltas)
    bootstrap_rows.append(
        {
            "variant": variant,
            "mean_MCC_delta": float(values.mean()),
            "ci025": float(np.quantile(values, 0.025)),
            "ci975": float(np.quantile(values, 0.975)),
            "probability_positive": float((values > 0).mean()),
        }
    )

pl.DataFrame(rows).write_csv(out / "metrics.csv")
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
pl.DataFrame(prediction_rows).write_parquet(out / "oof.parquet")
(out / "summary.json").write_text(json.dumps({"selection": cfg["selection"]}, indent=2) + "\n")
print(pl.DataFrame(rows), flush=True)
