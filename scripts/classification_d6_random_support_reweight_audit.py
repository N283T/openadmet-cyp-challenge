"""Post-stratify random-fold D6 OOF predictions to blind-test support."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(key: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[key])


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = rankdata(values[mask], method="average") / (mask.sum() + 1.0)
    return result


def weighted_calls(scores: np.ndarray, weights: np.ndarray) -> np.ndarray:
    order = np.argsort(-scores, kind="stable")
    target = cfg["call_rate"] * weights.sum()
    count = int(np.searchsorted(np.cumsum(weights[order]), target, side="left") + 1)
    calls = np.zeros(len(scores), dtype=bool)
    calls[order[:count]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray, weights: np.ndarray) -> float:
    positive = y == 1
    tpr = float(np.average(calls[positive], weights=weights[positive]))
    fpr = float(np.average(calls[~positive], weights=weights[~positive]))
    prevalence = cfg["target_prior"]
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


compounds = read("compounds").sort("compound_id")
labels = read("labels").filter(pl.col("endpoint") == cfg["endpoint"])
random_assignment = read("assignments").filter(pl.col("strategy") == "random").select(
    "compound_id", "fold"
)
data = labels.join(random_assignment, on="compound_id", validate="1:1").sort("compound_id")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
compound_ids = compounds["compound_id"].to_list()
index = {compound_id: row for row, compound_id in enumerate(compound_ids)}
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fingerprints = [
    generator.GetFingerprint(Chem.MolFromSmiles(smiles))
    for smiles in compounds["standardized_smiles"]
]
data_index = np.asarray([index[value] for value in data["compound_id"]])
test_index = np.asarray([index[value] for value in test["compound_id"]])
folds = data["fold"].to_numpy().astype(int)
y = data["value"].to_numpy().astype(int)

oof_support = np.empty(len(data), dtype=float)
for fold in np.unique(folds):
    train_fingerprints = [fingerprints[row] for row in data_index[folds != fold]]
    for row in np.flatnonzero(folds == fold):
        oof_support[row] = max(
            DataStructs.BulkTanimotoSimilarity(fingerprints[data_index[row]], train_fingerprints)
        )
full_train_fingerprints = [fingerprints[row] for row in data_index]
test_support = np.asarray(
    [
        max(DataStructs.BulkTanimotoSimilarity(fingerprints[row], full_train_fingerprints))
        for row in test_index
    ]
)

edges = np.asarray(cfg["support_bins"])
oof_bin = np.clip(np.digitize(oof_support, edges[1:-1]), 0, len(edges) - 2)
test_bin = np.clip(np.digitize(test_support, edges[1:-1]), 0, len(edges) - 2)
oof_counts = np.bincount(oof_bin, minlength=len(edges) - 1)
test_counts = np.bincount(test_bin, minlength=len(edges) - 1)
oof_share = oof_counts / oof_counts.sum()
test_share = test_counts / test_counts.sum()
if np.any((test_share > 0) & (oof_share == 0)):
    raise ValueError("Blind-test support mass exists in an empty random-OOF bin")
bin_weights = np.divide(test_share, oof_share, out=np.zeros_like(test_share), where=oof_share > 0)
weights = bin_weights[oof_bin]
weights /= weights.mean()

models: list[tuple[str, np.ndarray, str]] = []
for source, prefix in [("tabpfn_oof", "tabpfn35fast"), ("tree_oof", "public_tree")]:
    table = read(source).filter(pl.col("split_id") == "random")
    for method in sorted(table["method"].unique().to_list()):
        selected = table.filter(pl.col("method") == method).select(
            "compound_id", "fold", "y_true", "prediction"
        )
        joined = data.select("compound_id").join(
            selected, on="compound_id", validate="1:1", maintain_order="left"
        )
        if not np.array_equal(joined["fold"].to_numpy(), folds):
            raise ValueError(f"Fold mismatch for {prefix}::{method}")
        models.append((f"{prefix}::{method}", joined["prediction"].to_numpy(), "random"))

deployed = (
    read("deployed_oof")
    .filter(pl.col("model") == cfg["deployed_model"])
    .select("compound_id", "prediction")
)
joined = data.select("compound_id").join(
    deployed, on="compound_id", validate="1:1", maintain_order="left"
)
models.append((cfg["deployed_model"], joined["prediction"].to_numpy(), "umap_reference"))

metric_rows = []
bin_rows = []
for name, raw, training_split in models:
    scores = fold_rank(raw, folds)
    calls = weighted_calls(scores, weights)
    metric_rows.append(
        {
            "model": name,
            "training_split": training_split,
            "weighted_AUROC": float(roc_auc_score(y, scores, sample_weight=weights)),
            "weighted_AP": float(average_precision_score(y, scores, sample_weight=weights)),
            "weighted_MCC_target_prior": expected_mcc(y, calls, weights),
            "selected_rows": int(calls.sum()),
            "selected_weight": float(weights[calls].sum()),
        }
    )
    for bin_index in range(len(edges) - 1):
        mask = oof_bin == bin_index
        positives = y[mask] == 1
        bin_rows.append(
            {
                "model": name,
                "support_low": float(edges[bin_index]),
                "support_high": float(edges[bin_index + 1]),
                "n": int(mask.sum()),
                "positives": int(positives.sum()),
                "called": int(calls[mask].sum()),
                "called_positive": int((calls[mask] & positives).sum()),
            }
        )

support_rows = []
for bin_index in range(len(edges) - 1):
    support_rows.append(
        {
            "support_low": float(edges[bin_index]),
            "support_high": float(edges[bin_index + 1]),
            "oof_count": int(oof_counts[bin_index]),
            "oof_share": float(oof_share[bin_index]),
            "test_count": int(test_counts[bin_index]),
            "test_share": float(test_share[bin_index]),
            "importance_weight": float(bin_weights[bin_index]),
        }
    )

metrics = pl.DataFrame(metric_rows).sort("weighted_MCC_target_prior", descending=True)
metrics.write_csv(out / "weighted-metrics.csv")
pl.DataFrame(support_rows).write_csv(out / "support-histogram.csv")
pl.DataFrame(bin_rows).write_csv(out / "selected-by-support.csv")
data.select("compound_id", "fold").with_columns(
    pl.Series("support", oof_support),
    pl.Series("support_bin", oof_bin),
    pl.Series("importance_weight", weights),
).write_parquet(out / "random-oof-support.parquet")
test.with_columns(
    pl.Series("support", test_support), pl.Series("support_bin", test_bin)
).write_parquet(out / "test-support.parquet")
(out / "summary.json").write_text(json.dumps({"protocol": cfg["protocol"]}, indent=2) + "\n")
print(metrics.head(20), flush=True)
print(pl.DataFrame(support_rows), flush=True)
