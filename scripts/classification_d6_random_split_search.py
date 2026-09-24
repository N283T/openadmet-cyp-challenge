"""Find D6 random folds whose molecular support is closest to blind test."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import wasserstein_distance
from sklearn.model_selection import KFold, StratifiedKFold

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(key: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[key])


compounds = read("compounds").sort("compound_id")
labels = (
    read("labels")
    .filter(pl.col("endpoint") == cfg["endpoint"])
    .select("compound_id", pl.col("value").cast(pl.Int64).alias("y"))
    .sort("compound_id")
)
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
compound_ids = compounds["compound_id"].to_list()
index = {compound_id: row for row, compound_id in enumerate(compound_ids)}
generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=cfg["fingerprint_radius"], fpSize=cfg["fingerprint_size"]
)
fingerprints = [
    generator.GetFingerprint(Chem.MolFromSmiles(smiles))
    for smiles in compounds["standardized_smiles"]
]
label_index = np.asarray([index[value] for value in labels["compound_id"]])
test_index = np.asarray([index[value] for value in test["compound_id"]])
label_fingerprints = [fingerprints[row] for row in label_index]

n = len(labels)
similarities = np.empty((n, n), dtype=np.float32)
for row, fingerprint in enumerate(label_fingerprints):
    similarities[row] = DataStructs.BulkTanimotoSimilarity(fingerprint, label_fingerprints)
np.fill_diagonal(similarities, -1.0)
neighbor_order = np.argsort(-similarities, axis=1, kind="stable").astype(np.int32)
ordered_similarities = np.take_along_axis(similarities, neighbor_order, axis=1)
test_support = np.asarray(
    [
        max(DataStructs.BulkTanimotoSimilarity(fingerprints[row], label_fingerprints))
        for row in test_index
    ]
)


def make_folds(seed: int) -> np.ndarray:
    folds = np.full(n, -1, dtype=int)
    if cfg.get("splitter", "kfold") == "stratified_kfold":
        iterator = StratifiedKFold(
            cfg["n_folds"], shuffle=True, random_state=seed
        ).split(np.arange(n), labels["y"].to_numpy())
    else:
        iterator = KFold(cfg["n_folds"], shuffle=True, random_state=seed).split(np.arange(n))
    for fold, (_, validation) in enumerate(iterator):
        folds[validation] = fold
    return folds


def fold_support(folds: np.ndarray) -> np.ndarray:
    other_fold = folds[neighbor_order] != folds[:, None]
    first = np.argmax(other_fold, axis=1)
    if not other_fold[np.arange(n), first].all():
        raise ValueError("No cross-fold neighbor found")
    return ordered_similarities[np.arange(n), first]


search_rows = []
folds_by_seed = {}
support_by_seed = {}
for seed in range(cfg["candidate_seeds"]):
    folds = make_folds(seed)
    support = fold_support(folds)
    folds_by_seed[seed] = folds
    support_by_seed[seed] = support
    search_rows.append(
        {
            "seed": seed,
            "wasserstein_to_test": float(wasserstein_distance(support, test_support)),
            "median": float(np.median(support)),
            "mean": float(support.mean()),
            "fraction_ge_0_5": float((support >= 0.5).mean()),
            "fraction_ge_0_6": float((support >= 0.6).mean()),
        }
    )

search = pl.DataFrame(search_rows).sort(["wasserstein_to_test", "seed"])
selected = search.head(cfg["selected_seeds"])["seed"].to_list()
if cfg["baseline_seed"] not in selected:
    selected.append(cfg["baseline_seed"])

y = labels["y"].to_numpy()
assignment_rows = []
balance_rows = []
support_rows = []
for seed in selected:
    folds = folds_by_seed[seed]
    support = support_by_seed[seed]
    strategy = f"random_seed_{seed}"
    for row, compound_id in enumerate(labels["compound_id"]):
        assignment_rows.append(
            {
                "strategy": strategy,
                "compound_id": compound_id,
                "fold": int(folds[row]),
                "group": str(row),
            }
        )
        support_rows.append(
            {
                "strategy": strategy,
                "compound_id": compound_id,
                "fold": int(folds[row]),
                "support": float(support[row]),
            }
        )
    for fold in range(cfg["n_folds"]):
        mask = folds == fold
        balance_rows.append(
            {
                "strategy": strategy,
                "fold": fold,
                "n": int(mask.sum()),
                "positives": int(y[mask].sum()),
                "positive_rate": float(y[mask].mean()),
            }
        )

search.write_csv(out / "seed-search.csv")
pl.DataFrame(assignment_rows).write_parquet(out / "assignments.parquet")
pl.DataFrame(balance_rows).write_csv(out / "fold-balance.csv")
pl.DataFrame(support_rows).write_parquet(out / "oof-support.parquet")
test.with_columns(pl.Series("support", test_support)).write_parquet(out / "test-support.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "selection": cfg["selection"],
            "selected_seeds": selected,
            "test_support_median": float(np.median(test_support)),
            "test_support_fraction_ge_0_6": float((test_support >= 0.6).mean()),
        },
        indent=2,
    )
    + "\n"
)
print(search.head(10), flush=True)
print(pl.DataFrame(balance_rows), flush=True)
