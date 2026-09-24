"""Audit UMAP-prespecified new-family D6 subsets on scaffold folds."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

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


def fixed_calls(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    calls = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(rows))))
        calls[rows[np.argsort(-values[rows], kind="stable")[:count]]] = True
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


sources = {
    "lgbm": read("lgbm_oof"),
    "readout": read("readout_oof"),
    "tabpfn": read("tabpfn_oof"),
}
reference_spec = next(member for member in cfg["members"] if member["name"] == cfg["reference"])
reference_table = (
    sources[reference_spec["source"]]
    .filter(
        (pl.col("endpoint") == cfg["endpoint"])
        & (pl.col("split_id") == cfg["split_id"])
        & (pl.col("method") == reference_spec["method"])
    )
    .sort("compound_id")
)
ids = reference_table.select("compound_id")
folds = reference_table["fold"].to_numpy().astype(int)
y = reference_table["y_true"].to_numpy().astype(int)
groups = reference_table["group"].cast(pl.String).to_numpy()


def align(spec: dict) -> np.ndarray:
    table = sources[spec["source"]].filter(
        (pl.col("endpoint") == cfg["endpoint"])
        & (pl.col("split_id") == cfg["split_id"])
        & (pl.col("method") == spec["method"])
    )
    joined = ids.join(
        table.select("compound_id", "prediction"),
        on="compound_id",
        validate="1:1",
        maintain_order="left",
    )
    if joined["prediction"].null_count() or not joined["prediction"].is_finite().all():
        raise ValueError(f"Alignment failed for {spec['name']}")
    return fold_rank(joined["prediction"].to_numpy(), folds)


names = [member["name"] for member in cfg["members"]]
x = np.column_stack([align(member) for member in cfg["members"]])
index = {name: position for position, name in enumerate(names)}
predictions = {"family_control": x[:, index[cfg["reference"]]]}
for subset_name, subset in cfg["subsets"].items():
    members = list(range(len(names))) if subset == "all" else [index[name] for name in subset]
    predictions[subset_name] = x[:, members].mean(axis=1)

metric_rows, calls = [], {}
for method, values in predictions.items():
    selected = fixed_calls(values, folds)
    calls[method] = selected
    metric_rows.append(
        {
            "method": method,
            "AUROC": float(roc_auc_score(y, values)),
            "AP": float(average_precision_score(y, values)),
            "MCC_observed": float(matthews_corrcoef(y, selected)),
            "MCC_target_prior": expected_mcc(y, selected),
            "positive_calls": int(selected.sum()),
            "swaps_from_control": int(
                np.count_nonzero(selected != fixed_calls(predictions["family_control"], folds))
            ),
        }
    )

rng = np.random.default_rng(cfg["bootstrap_seed"])
unique_groups = np.unique(groups)
group_rows = {group: np.flatnonzero(groups == group) for group in unique_groups}
bootstrap_rows = []
for candidate in cfg["subsets"]:
    delta = []
    for _ in range(cfg["bootstrap_replicates"]):
        sampled = rng.choice(unique_groups, len(unique_groups), replace=True)
        rows = np.concatenate([group_rows[group] for group in sampled])
        delta.append(
            expected_mcc(y[rows], calls[candidate][rows])
            - expected_mcc(y[rows], calls["family_control"][rows])
        )
    values = np.asarray(delta)
    bootstrap_rows.append(
        {
            "candidate": candidate,
            "control": "family_control",
            "mean_delta": float(values.mean()),
            "ci025": float(np.quantile(values, 0.025)),
            "ci975": float(np.quantile(values, 0.975)),
            "probability_positive": float((values > 0).mean()),
        }
    )

metrics = pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "protocol": cfg["protocol"],
            "candidate_created": False,
            "uploaded": False,
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
