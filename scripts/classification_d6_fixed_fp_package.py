"""Audit and package the fold-selected fixed-Morgan classical D6 candidate."""

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
rate = 49 / 750


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


def calls(values, folds):
    binary = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = int(round(rate * len(positions)))
        order = np.argsort(-values[positions], kind="stable")
        binary[positions[order[:count]]] = True
    return binary


oof = read("oof").sort("compound_id")
selection = pl.read_csv(inp / cfg["selection"])
top1 = selection.filter(pl.col("variant") == "top1").row(0, named=True)
fold_methods = json.loads(top1["fold_methods"])
live = oof.filter(pl.col("method") == "live")
y = live["y_true"].to_numpy().astype(int)
folds = live["fold"].to_numpy()
groups = live["group"].to_numpy()
crossfit_prediction = np.empty(len(y))
for fold, selected in enumerate(fold_methods):
    assert len(selected) == 1
    table = oof.filter(pl.col("method") == selected[0])
    assert table["compound_id"].equals(live["compound_id"])
    mask = folds == fold
    crossfit_prediction[mask] = table["prediction"].to_numpy()[mask]
live_calls = calls(live["prediction"].to_numpy(), folds)
crossfit_calls = calls(crossfit_prediction, folds)

fold_rows = []
for fold in [-1, 0, 1, 2, 3, 4]:
    mask = np.ones(len(y), dtype=bool) if fold == -1 else folds == fold
    fold_rows.append(
        {
            "fold": fold,
            "selected": "crossfit" if fold == -1 else fold_methods[fold][0],
            "live_expected_MCC": expected_mcc(y[mask], live_calls[mask], cfg["target_prior"]),
            "candidate_expected_MCC": expected_mcc(
                y[mask], crossfit_calls[mask], cfg["target_prior"]
            ),
        }
    )
pl.DataFrame(fold_rows).write_csv(out / "fold-comparison.csv")

group_indices = [np.flatnonzero(groups == group) for group in np.unique(groups)]
rng = np.random.default_rng(20260913)
differences = []
for _ in range(10000):
    selected = rng.integers(0, len(group_indices), len(group_indices))
    indices = np.concatenate([group_indices[position] for position in selected])
    differences.append(
        expected_mcc(y[indices], crossfit_calls[indices], cfg["target_prior"])
        - expected_mcc(y[indices], live_calls[indices], cfg["target_prior"])
    )
bootstrap = {
    "difference": expected_mcc(y, crossfit_calls, cfg["target_prior"])
    - expected_mcc(y, live_calls, cfg["target_prior"]),
    "low": float(np.quantile(differences, 0.025)),
    "high": float(np.quantile(differences, 0.975)),
    "resamples": 10000,
}

deployment = read("deployment").filter(pl.col("method") == "top1").sort("compound_id")
assert deployment.height == 750 and int(deployment["binary"].sum()) == 49
aliases = (
    read("aliases")
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
live_submission = read("live_submission").sort("Molecule_Name")
live_with_ids = live_submission.join(aliases, on="Molecule_Name", validate="1:1")
live_d6 = live_with_ids.sort("compound_id")[endpoint].to_numpy().astype(bool)
candidate_d6 = deployment["binary"].to_numpy()
overlap = {
    "both_positive": int((live_d6 & candidate_d6).sum()),
    "live_only": int((live_d6 & ~candidate_d6).sum()),
    "candidate_only": int((~live_d6 & candidate_d6).sum()),
    "both_negative": int((~live_d6 & ~candidate_d6).sum()),
}
candidate = (
    live_with_ids.drop(endpoint)
    .join(
        deployment.select("compound_id", pl.col("binary").alias(endpoint)),
        on="compound_id",
        validate="1:1",
    )
    .select("SMILES", "Molecule_Name", *ENDPOINTS["classification"])
    .sort("Molecule_Name")
)
candidate.write_parquet(out / "fixed-morgan-model-selection.parquet")
candidate.write_csv(out / "fixed-morgan-model-selection.csv")
report = validate(
    out / "fixed-morgan-model-selection.parquet", "classification", inp / cfg["reference"]
)
(out / "summary.json").write_text(
    json.dumps(
        {
            "bootstrap": bootstrap,
            "test_overlap": overlap,
            "deployment_method": json.loads(top1["deployment_methods"]),
            "fold_methods": fold_methods,
            "positive_calls": 49,
            "validation": report,
            "uploaded": False,
        },
        indent=2,
    )
    + "\n"
)
print((out / "summary.json").read_text(), flush=True)
