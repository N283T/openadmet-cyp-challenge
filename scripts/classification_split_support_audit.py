"""Compare candidate CV splits with blind-test structural support."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import wasserstein_distance

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def describe(values: np.ndarray) -> dict[str, float]:
    return {
        "mean": float(values.mean()),
        "p10": float(np.quantile(values, 0.10)),
        "p25": float(np.quantile(values, 0.25)),
        "median": float(np.median(values)),
        "p75": float(np.quantile(values, 0.75)),
        "p90": float(np.quantile(values, 0.90)),
        "fraction_lt_0_3": float((values < 0.3).mean()),
        "fraction_lt_0_4": float((values < 0.4).mean()),
        "fraction_ge_0_6": float((values >= 0.6).mean()),
    }


def nearest_support(query: np.ndarray, pool: np.ndarray, fingerprints: list) -> np.ndarray:
    pool_fingerprints = [fingerprints[index] for index in pool]
    result = np.empty(len(query), dtype=float)
    for position, index in enumerate(query):
        result[position] = max(
            DataStructs.BulkTanimotoSimilarity(fingerprints[index], pool_fingerprints)
        )
    return result


compounds = read("compounds").sort("compound_id")
labels = read("labels")
test_ids = set(read("test_mapping")["compound_id"].to_list())
assignments = read("assignments")
compound_ids = compounds["compound_id"].to_list()
id_to_index = {compound_id: row for row, compound_id in enumerate(compound_ids)}
test_index = np.asarray([id_to_index[value] for value in compound_ids if value in test_ids])

generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=cfg["fingerprint"]["radius"], fpSize=cfg["fingerprint"]["size"]
)
fingerprints = [
    generator.GetFingerprint(Chem.MolFromSmiles(value))
    for value in compounds["standardized_smiles"]
]

summary_rows = []
distribution_rows = []
fold_rows = []
for endpoint in cfg["endpoints"]:
    endpoint_labels = labels.filter(pl.col("endpoint") == endpoint).select("compound_id", "value")
    endpoint_ids = endpoint_labels["compound_id"].to_list()
    endpoint_index = np.asarray([id_to_index[value] for value in endpoint_ids])
    target_support = nearest_support(test_index, endpoint_index, fingerprints)
    target_desc = describe(target_support)
    for value in target_support:
        distribution_rows.append(
            {
                "endpoint": endpoint,
                "strategy": "blind_test_target",
                "fold": -1,
                "support": float(value),
            }
        )

    leave_one_out = np.empty(len(endpoint_index), dtype=float)
    endpoint_fingerprints = [fingerprints[index] for index in endpoint_index]
    for row, fingerprint in enumerate(endpoint_fingerprints):
        similarities = DataStructs.BulkTanimotoSimilarity(fingerprint, endpoint_fingerprints)
        leave_one_out[row] = max(
            value for position, value in enumerate(similarities) if position != row
        )
    ceiling_desc = describe(leave_one_out)
    summary_rows.append(
        {
            "endpoint": endpoint,
            "strategy": "leave_one_out_ceiling",
            "n_oof": len(leave_one_out),
            "n_test": len(target_support),
            "wasserstein_to_test": float(wasserstein_distance(leave_one_out, target_support)),
            "absolute_median_gap": abs(ceiling_desc["median"] - target_desc["median"]),
            "absolute_lt_0_4_gap": abs(
                ceiling_desc["fraction_lt_0_4"] - target_desc["fraction_lt_0_4"]
            ),
            **{f"oof_{key}": value for key, value in ceiling_desc.items()},
            **{f"test_{key}": value for key, value in target_desc.items()},
        }
    )
    for value in leave_one_out:
        distribution_rows.append(
            {
                "endpoint": endpoint,
                "strategy": "leave_one_out_ceiling",
                "fold": -1,
                "support": float(value),
            }
        )

    for strategy in cfg["strategies"]:
        split = assignments.filter(pl.col("strategy") == strategy).select(
            "compound_id", "fold", "group"
        )
        data = endpoint_labels.join(split, on="compound_id", validate="1:1").sort("compound_id")
        indices = np.asarray([id_to_index[value] for value in data["compound_id"]])
        folds = data["fold"].to_numpy()
        values = data["value"].to_numpy().astype(int)
        support = np.empty(len(data), dtype=float)
        for fold in sorted(set(folds)):
            heldout = np.flatnonzero(folds == fold)
            training = np.flatnonzero(folds != fold)
            support[heldout] = nearest_support(indices[heldout], indices[training], fingerprints)
            fold_rows.append(
                {
                    "endpoint": endpoint,
                    "strategy": strategy,
                    "fold": int(fold),
                    "n": int(len(heldout)),
                    "positives": int(values[heldout].sum()),
                    "positive_rate": float(values[heldout].mean()),
                    "groups": int(data.filter(pl.col("fold") == fold)["group"].n_unique()),
                    **describe(support[heldout]),
                }
            )
        for fold, value in zip(folds, support, strict=True):
            distribution_rows.append(
                {
                    "endpoint": endpoint,
                    "strategy": strategy,
                    "fold": int(fold),
                    "support": float(value),
                }
            )
        desc = describe(support)
        summary_rows.append(
            {
                "endpoint": endpoint,
                "strategy": strategy,
                "n_oof": len(support),
                "n_test": len(target_support),
                "wasserstein_to_test": float(wasserstein_distance(support, target_support)),
                "absolute_median_gap": abs(desc["median"] - target_desc["median"]),
                "absolute_lt_0_4_gap": abs(
                    desc["fraction_lt_0_4"] - target_desc["fraction_lt_0_4"]
                ),
                **{f"oof_{key}": value for key, value in desc.items()},
                **{f"test_{key}": value for key, value in target_desc.items()},
            }
        )

summary = pl.DataFrame(summary_rows).sort("endpoint", "wasserstein_to_test")
summary.write_csv(out / "summary.csv")
pl.DataFrame(fold_rows).sort("endpoint", "strategy", "fold").write_csv(out / "fold-summary.csv")
pl.DataFrame(distribution_rows).write_parquet(out / "support-distributions.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "fingerprint": cfg["fingerprint"],
            "test_labels_used": False,
            "selection_rule": (
                "Choose a prespecified grouped secondary strategy by the smallest "
                "endpoint-averaged Wasserstein distance to blind-test support, subject "
                "to nondegenerate endpoint label counts in every fold."
            ),
            "strategies": cfg["strategies"],
        },
        indent=2,
    )
    + "\n"
)
print(summary.write_csv(), flush=True)
