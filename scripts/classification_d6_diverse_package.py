"""Audit and package diverse CYP2D6 fixed-rate candidates with two CYP3A4 choices."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.tracks import ENDPOINTS, validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoint = "CYP2D6_is_TDI"
target_prior = cfg["target_prior"]


def read(name):
    return pl.read_parquet(inp / cfg[name])


def expected_mcc(y, binary, prevalence):
    positive = y == 1
    tpr = float(binary[positive].mean())
    fpr = float(binary[~positive].mean())
    tp = prevalence * tpr
    fn = prevalence * (1 - tpr)
    fp = (1 - prevalence) * fpr
    tn = (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator > 0 else 0.0


def calls_by_fold(values, folds, rate):
    calls = np.zeros(len(values), dtype=bool)
    for fold in range(5):
        positions = np.flatnonzero(folds == fold)
        count = int(round(rate * len(positions)))
        order = np.argsort(-values[positions], kind="stable")
        calls[positions[order[:count]]] = True
    return calls


def top_calls(values, count):
    calls = np.zeros(len(values), dtype=bool)
    calls[np.argsort(-values, kind="stable")[:count]] = True
    return calls


oof = read("oof").sort("compound_id")
test_scores = read("test_scores").sort("compound_id")
reference_oof = oof.filter(pl.col("method") == "current")
assert reference_oof.height == 1497
y = reference_oof["y_true"].to_numpy().astype(int)
folds = reference_oof["fold"].to_numpy()
groups = reference_oof["group"].to_numpy()

candidate_specs = cfg["candidates"]
evaluation_rows = []
binary_oof = {}
for name, spec in candidate_specs.items():
    table = oof.filter(pl.col("method") == spec["method"])
    assert table["compound_id"].equals(reference_oof["compound_id"])
    calls = calls_by_fold(table["prediction"].to_numpy(), folds, spec["positives"] / 750)
    binary_oof[name] = calls
    for fold in [-1, 0, 1, 2, 3, 4]:
        mask = np.ones(len(y), dtype=bool) if fold == -1 else folds == fold
        evaluation_rows.append(
            {
                "candidate": name,
                "fold": fold,
                "n": int(mask.sum()),
                "positives": int(calls[mask].sum()),
                "expected_target_MCC": expected_mcc(y[mask], calls[mask], target_prior),
            }
        )
pl.DataFrame(evaluation_rows).write_csv(out / "fold-metrics.csv")

group_indices = [np.flatnonzero(groups == group) for group in np.unique(groups)]
rng = np.random.default_rng(20260913)
bootstrap_rows = []
reference_calls = binary_oof["current_43"]
difference_distributions = {
    name: [] for name in binary_oof if name != "current_43"
}
for _ in range(10000):
    selected = rng.integers(0, len(group_indices), len(group_indices))
    indices = np.concatenate([group_indices[position] for position in selected])
    reference_score = expected_mcc(y[indices], reference_calls[indices], target_prior)
    for name in difference_distributions:
        difference_distributions[name].append(
            expected_mcc(y[indices], binary_oof[name][indices], target_prior)
            - reference_score
        )
for name, calls in binary_oof.items():
    if name == "current_43":
        continue
    differences = difference_distributions[name]
    bootstrap_rows.append(
        {
            "candidate": name,
            "reference": "current_43",
            "difference": expected_mcc(y, calls, target_prior)
            - expected_mcc(y, reference_calls, target_prior),
            "low": float(np.quantile(differences, 0.025)),
            "high": float(np.quantile(differences, 0.975)),
            "resamples": 10000,
        }
    )
pl.DataFrame(bootstrap_rows).write_csv(out / "conditional-bootstrap.csv")

aliases = (
    read("aliases")
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
reference = pl.read_csv(inp / cfg["reference"])
identities = reference.select("SMILES", "Molecule_Name").join(
    aliases, on="Molecule_Name", validate="1:1"
)
assert identities.height == 750
a4_sources = {
    "caruana_a4": read("caruana_a4"),
    "old_a4": read("old_a4"),
}
reports = []
overlaps = []
current_test = test_scores.filter(pl.col("method") == "current")
current_test_calls = top_calls(current_test["prediction"].to_numpy(), 43)
for name, spec in candidate_specs.items():
    if name == "current_43":
        continue
    scores = test_scores.filter(pl.col("method") == spec["method"])
    assert scores["compound_id"].equals(current_test["compound_id"])
    calls = top_calls(scores["prediction"].to_numpy(), spec["positives"])
    overlaps.append(
        {
            "candidate": name,
            "both_positive": int((calls & current_test_calls).sum()),
            "current_only": int((~calls & current_test_calls).sum()),
            "candidate_only": int((calls & ~current_test_calls).sum()),
            "both_negative": int((~calls & ~current_test_calls).sum()),
        }
    )
    d6 = scores.select("compound_id").with_columns(pl.Series(endpoint, calls))
    for a4_name, a4_source in a4_sources.items():
        candidate = (
            identities.join(d6, on="compound_id", validate="1:1")
            .join(
                a4_source.select("Molecule_Name", "CYP3A4_is_TDI"),
                on="Molecule_Name",
                validate="1:1",
            )
            .select("SMILES", "Molecule_Name", *ENDPOINTS["classification"])
            .sort("Molecule_Name")
        )
        candidate_name = f"{name}_{a4_name}"
        candidate.write_parquet(out / f"{candidate_name}.parquet")
        candidate.write_csv(out / f"{candidate_name}.csv")
        reports.append(
            {
                "candidate": candidate_name,
                "D6_positives": int(candidate[endpoint].sum()),
                "A4_positives": int(candidate["CYP3A4_is_TDI"].sum()),
                "validation": validate(
                    out / f"{candidate_name}.parquet", "classification", inp / cfg["reference"]
                ),
            }
        )
pl.DataFrame(overlaps).write_csv(out / "test-overlap.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "candidates": reports,
            "comparison": "Conditional screen after choosing candidates on these OOF results",
            "bootstrap": "UMAP-group resampling; fixed original OOF calls",
            "uploaded": False,
        },
        indent=2,
    )
    + "\n"
)
print(json.dumps(reports, indent=2), flush=True)
