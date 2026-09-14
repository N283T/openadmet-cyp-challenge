"""Gate the independent FDA TDI block by similarity to its source domain."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import rankdata
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


def fixed_rate_calls(scores: np.ndarray, folds: np.ndarray) -> np.ndarray:
    calls = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(positions))))
        calls[positions[np.argsort(-scores[positions], kind="stable")[:count]]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray) -> float:
    positive = y == 1
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    prevalence = cfg["target_prior"]
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


external = read("external_records").filter(
    pl.col("value").is_in([0, 1]) & pl.col("challenge_compound_id").is_null()
)
compounds = read("compounds").sort("compound_id")
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
external_fp = [
    generator.GetFingerprint(Chem.MolFromSmiles(value)) for value in external["canonical_smiles"]
]
challenge_fp = [
    generator.GetFingerprint(Chem.MolFromSmiles(value))
    for value in compounds["standardized_smiles"]
]
max_similarity = np.asarray(
    [max(DataStructs.BulkTanimotoSimilarity(fp, external_fp)) for fp in challenge_fp]
)
id_to_index = {value: row for row, value in enumerate(compounds["compound_id"])}
prior_oof = read("prior_oof")
prior_test = read("prior_test")

metric_rows, oof_rows, test_rows = [], [], []
for strategy in cfg["strategies"]:
    base = prior_oof.filter(
        (pl.col("strategy") == strategy) & (pl.col("method") == "challenge_block")
    ).sort("compound_id")
    external_oof = prior_oof.filter(
        (pl.col("strategy") == strategy) & (pl.col("method") == "fda_lr")
    ).sort("compound_id")
    deploy = prior_test.filter(
        (pl.col("strategy") == strategy) & (pl.col("method") == "challenge_block")
    ).sort("compound_id")
    external_deploy = prior_test.filter(
        (pl.col("strategy") == strategy) & (pl.col("method") == "fda_lr")
    ).sort("compound_id")
    folds = base["fold"].to_numpy().astype(int)
    y = base["y_true"].to_numpy().astype(int)
    oof_index = np.asarray([id_to_index[value] for value in base["compound_id"]])
    test_index = np.asarray([id_to_index[value] for value in deploy["compound_id"]])
    base_rank = fold_rank(base["prediction"].to_numpy(), folds)
    external_rank = fold_rank(external_oof["prediction"].to_numpy(), folds)
    base_test_rank = unit_rank(deploy["prediction"].to_numpy())
    external_test_rank = unit_rank(external_deploy["prediction"].to_numpy())
    gates = {
        f"threshold_{threshold}": (max_similarity >= threshold).astype(float)
        for threshold in cfg["thresholds"]
    }
    gates["continuous"] = np.clip(
        (max_similarity - cfg["continuous_low"]) / (cfg["continuous_high"] - cfg["continuous_low"]),
        0,
        1,
    )
    candidates = {"challenge_block": (base_rank, base_test_rank)}
    for gate_name, gate in gates.items():
        for alpha in cfg["alphas"]:
            oof_gate = gate[oof_index]
            test_gate = gate[test_index]
            candidates[f"{gate_name}_alpha_{alpha}"] = (
                base_rank + alpha * oof_gate * (external_rank - base_rank),
                base_test_rank + alpha * test_gate * (external_test_rank - base_test_rank),
            )
    for method, (prediction, deploy_prediction) in candidates.items():
        calls = fixed_rate_calls(prediction, folds)
        metric_rows.append(
            {
                "endpoint": cfg["endpoint"],
                "strategy": strategy,
                "method": method,
                "AUROC": float(roc_auc_score(y, prediction)),
                "AP": float(average_precision_score(y, prediction)),
                "MCC_observed": float(matthews_corrcoef(y, calls)),
                "MCC_target_prior": expected_mcc(y, calls),
            }
        )
        for row, value in enumerate(prediction):
            oof_rows.append(
                {
                    "compound_id": base["compound_id"][row],
                    "endpoint": cfg["endpoint"],
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
                    "endpoint": cfg["endpoint"],
                    "strategy": strategy,
                    "method": method,
                    "prediction": float(value),
                }
            )

metrics = pl.DataFrame(metric_rows).sort("strategy", "MCC_target_prior", descending=[False, True])
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(oof_rows).write_parquet(out / "oof.parquet")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
support = compounds.select("compound_id").with_columns(
    pl.Series("fda_max_tanimoto", max_similarity)
)
support.write_parquet(out / "support.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "median_max_tanimoto": float(np.median(max_similarity)),
            "test_median_max_tanimoto": float(np.median(max_similarity[test_index])),
            "scope": "Prespecified similarity-gated reuse of the independent FDA LR block.",
            "selection": "No candidate packaging or submission.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
