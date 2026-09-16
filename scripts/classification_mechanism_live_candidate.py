"""Audit and package the fixed endpoint-specific live/mechanism TDI candidate."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

from cyp_challenge.tracks import ENDPOINTS, validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
D6 = "CYP2D6_is_TDI"
A4 = "CYP3A4_is_TDI"


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


def fixed_calls(values: np.ndarray, count: int, folds: np.ndarray | None = None) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    if folds is None:
        result[np.argsort(-values, kind="stable")[:count]] = True
        return result
    rate = count / 750
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        fold_count = int(round(rate * len(rows)))
        result[rows[np.argsort(-values[rows], kind="stable")[:fold_count]]] = True
    return result


def expected_mcc(y: np.ndarray, calls: np.ndarray, prevalence: float) -> float:
    positive = y == 1
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def aligned_long(
    table: pl.DataFrame, endpoint: str, method: str, reference: pl.DataFrame
) -> pl.DataFrame:
    filtered = table.filter(
        (pl.col("endpoint") == endpoint)
        & (pl.col("method") == method)
        & ((pl.col("strategy") == "umap") if "strategy" in table.columns else pl.lit(True))
    ).select("compound_id", "prediction")
    joined = reference.select("compound_id").join(filtered, on="compound_id", validate="1:1")
    if joined.height != reference.height:
        raise ValueError(f"Incomplete {endpoint}/{method} alignment")
    return joined


# Reconstruct the current D6 OOF cross-fit policy and use its deployed top1 score.
d6_bank = read("d6_current_oof").sort("compound_id")
d6_live = d6_bank.filter(pl.col("method") == "live").sort("compound_id")
d6_selection = pl.read_csv(inp / cfg["d6_selection"]).filter(pl.col("variant") == "top1").row(
    0, named=True
)
d6_fold_methods = json.loads(d6_selection["fold_methods"])
d6_folds = d6_live["fold"].to_numpy().astype(int)
d6_current_raw = np.empty(d6_live.height, dtype=float)
for fold, methods in enumerate(d6_fold_methods):
    if len(methods) != 1:
        raise ValueError("D6 top1 policy must select exactly one method per fold")
    selected = d6_bank.filter(pl.col("method") == methods[0]).sort("compound_id")
    if not selected.select("compound_id").equals(d6_live.select("compound_id")):
        raise ValueError("D6 OOF method keys differ")
    mask = d6_folds == fold
    d6_current_raw[mask] = selected["prediction"].to_numpy()[mask]
d6_current_rank = fold_rank(d6_current_raw, d6_folds)
d6_mechanism = aligned_long(
    read("d6_mechanism_oof"), D6, cfg["d6_mechanism_method"], d6_live
)
d6_mechanism_rank = fold_rank(d6_mechanism["prediction"].to_numpy(), d6_folds)

d6_test = (
    read("d6_current_test")
    .filter(pl.col("method") == "top1")
    .select("compound_id", "prediction")
    .sort("compound_id")
)
d6_mechanism_test = aligned_long(
    read("d6_mechanism_test"), D6, cfg["d6_mechanism_method"], d6_test
)
d6_current_test_rank = unit_rank(d6_test["prediction"].to_numpy())
d6_mechanism_test_rank = unit_rank(d6_mechanism_test["prediction"].to_numpy())

# Reconstruct the currently deployed A4 80/20 Caruana/PXR-Admetica OOF rank.
a4_caruana = read("a4_caruana_oof").sort("compound_id")
a4_folds = a4_caruana["fold"].to_numpy().astype(int)
a4_caruana_rank = fold_rank(a4_caruana["caruana_auc"].to_numpy(), a4_folds)
a4_pxr = aligned_long(read("a4_pxr_oof"), A4, cfg["a4_pxr_method"], a4_caruana)
a4_pxr_rank = fold_rank(a4_pxr["prediction"].to_numpy(), a4_folds)
a4_current_rank = (
    (1 - cfg["a4_pxr_weight"]) * a4_caruana_rank
    + cfg["a4_pxr_weight"] * a4_pxr_rank
)
a4_mechanism = aligned_long(
    read("a4_mechanism_oof"), A4, cfg["a4_mechanism_method"], a4_caruana
)
a4_group_source = read("a4_mechanism_oof").filter(
    (pl.col("endpoint") == A4)
    & (pl.col("method") == cfg["a4_mechanism_method"])
    & (pl.col("strategy") == "umap")
).select("compound_id", "group")
a4_groups = (
    a4_caruana.select("compound_id")
    .join(a4_group_source, on="compound_id", validate="1:1")
    ["group"]
    .cast(pl.String)
    .to_numpy()
)
a4_mechanism_rank = fold_rank(a4_mechanism["prediction"].to_numpy(), a4_folds)

a4_test = (
    read("a4_current_test")
    .filter(pl.col("method") == cfg["a4_current_test_method"])
    .select("compound_id", "prediction")
    .sort("compound_id")
)
a4_mechanism_test = aligned_long(
    read("a4_mechanism_test"), A4, cfg["a4_mechanism_method"], a4_test
)
a4_current_test_rank = unit_rank(a4_test["prediction"].to_numpy())
a4_mechanism_test_rank = unit_rank(a4_mechanism_test["prediction"].to_numpy())

weight = cfg["blend_weight"]
endpoint_data = {
    D6: {
        "ids": d6_live["compound_id"],
        "y": d6_live["y_true"].to_numpy().astype(int),
        "folds": d6_folds,
        "groups": d6_live["group"].cast(pl.String).to_numpy(),
        "current": d6_current_rank,
        "mechanism": d6_mechanism_rank,
        "test_ids": d6_test["compound_id"],
        "current_test": d6_current_test_rank,
        "mechanism_test": d6_mechanism_test_rank,
        "count": cfg["d6_positive_calls"],
    },
    A4: {
        "ids": a4_caruana["compound_id"],
        "y": a4_caruana["y_true"].to_numpy().astype(int),
        "folds": a4_folds,
        "groups": a4_groups,
        "current": a4_current_rank,
        "mechanism": a4_mechanism_rank,
        "test_ids": a4_test["compound_id"],
        "current_test": a4_current_test_rank,
        "mechanism_test": a4_mechanism_test_rank,
        "count": cfg["a4_positive_calls"],
    },
}

metric_rows = []
bootstrap_rows = []
test_score_rows = []
selected_test = {}
rng = np.random.default_rng(cfg["bootstrap_seed"])
for endpoint, data in endpoint_data.items():
    blend = (1 - weight) * data["current"] + weight * data["mechanism"]
    test_blend = (1 - weight) * data["current_test"] + weight * data["mechanism_test"]
    current_calls = fixed_calls(data["current"], data["count"], data["folds"])
    blend_calls = fixed_calls(blend, data["count"], data["folds"])
    current_test_calls = fixed_calls(data["current_test"], data["count"])
    blend_test_calls = fixed_calls(test_blend, data["count"])
    selected_test[endpoint] = pl.DataFrame(
        {"compound_id": data["test_ids"], endpoint: blend_test_calls}
    )
    prevalence = cfg["target_prior"][endpoint]
    for method, scores, calls in [
        ("current", data["current"], current_calls),
        (
            "mechanism",
            data["mechanism"],
            fixed_calls(data["mechanism"], data["count"], data["folds"]),
        ),
        ("current75_mechanism25", blend, blend_calls),
    ]:
        metric_rows.append(
            {
                "endpoint": endpoint,
                "method": method,
                "AUROC": float(roc_auc_score(data["y"], scores)),
                "AP": float(average_precision_score(data["y"], scores)),
                "MCC_observed": float(matthews_corrcoef(data["y"], calls)),
                "MCC_target_prior": expected_mcc(data["y"], calls, prevalence),
                "calls": int(calls.sum()),
            }
        )
    group_values = np.unique(data["groups"])
    group_rows = {group: np.flatnonzero(data["groups"] == group) for group in group_values}
    deltas = []
    for _ in range(cfg["bootstrap_replicates"]):
        sampled = rng.choice(group_values, size=len(group_values), replace=True)
        rows = np.concatenate([group_rows[group] for group in sampled])
        deltas.append(
            expected_mcc(data["y"][rows], blend_calls[rows], prevalence)
            - expected_mcc(data["y"][rows], current_calls[rows], prevalence)
        )
    bootstrap_rows.append(
        {
            "endpoint": endpoint,
            "delta_target_prior_MCC": expected_mcc(data["y"], blend_calls, prevalence)
            - expected_mcc(data["y"], current_calls, prevalence),
            "p025": float(np.quantile(deltas, 0.025)),
            "p975": float(np.quantile(deltas, 0.975)),
            "improvement_probability": float(np.mean(np.asarray(deltas) > 0)),
            "oof_spearman_current_mechanism": float(
                spearmanr(data["current"], data["mechanism"]).statistic
            ),
            "test_spearman_current_mechanism": float(
                spearmanr(data["current_test"], data["mechanism_test"]).statistic
            ),
            "test_current_only": int((current_test_calls & ~blend_test_calls).sum()),
            "test_candidate_only": int((~current_test_calls & blend_test_calls).sum()),
        }
    )
    for compound_id, score, call in zip(
        data["test_ids"], test_blend, blend_test_calls, strict=True
    ):
        test_score_rows.append(
            {
                "compound_id": compound_id,
                "endpoint": endpoint,
                "method": "current75_mechanism25",
                "prediction": float(score),
                "call": bool(call),
                "fold": -1,
            }
        )

aliases = (
    read("aliases")
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
candidate = read("current_candidate").sort("Molecule_Name")
for endpoint in [D6, A4]:
    replacement = aliases.join(selected_test[endpoint], on="compound_id", validate="1:1").select(
        "Molecule_Name", endpoint
    )
    candidate = candidate.drop(endpoint).join(replacement, on="Molecule_Name", validate="1:1")
candidate = candidate.select("SMILES", "Molecule_Name", *ENDPOINTS["classification"]).sort(
    "Molecule_Name"
)

candidate.write_parquet(out / "predictions.parquet")
candidate.write_csv(out / "predictions.csv")
pl.DataFrame(metric_rows).sort("endpoint", "method").write_csv(out / "metrics.csv")
pl.DataFrame(bootstrap_rows).sort("endpoint").write_csv(out / "conditional-bootstrap.csv")
pl.DataFrame(test_score_rows).write_parquet(out / "test-scores.parquet")
summary = {
    "selection": cfg["selection"],
    "blend_weight": weight,
    "positive_calls": {endpoint: int(candidate[endpoint].sum()) for endpoint in [D6, A4]},
    "validation": validate(out / "predictions.parquet", "classification", inp / cfg["reference"]),
    "uploaded": False,
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(pl.DataFrame(metric_rows).sort("endpoint", "method"), flush=True)
print(pl.DataFrame(bootstrap_rows).sort("endpoint"), flush=True)
print(json.dumps(summary, indent=2), flush=True)
