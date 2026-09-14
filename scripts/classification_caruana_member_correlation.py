"""Audit dependence and standalone performance of deployed A4 Caruana members."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
selections = json.loads((inp / cfg["selections"]).read_text())
selection = next(
    row
    for row in selections
    if row["endpoint"] == cfg["endpoint"] and row["recipe"] == cfg["recipe"] and row["fold"] == -1
)
weights = {name: value for name, value in selection["weights"].items() if value > 0}


def prediction_column(table: pl.DataFrame) -> str:
    return "probability" if "probability" in table.columns else "prediction"


def read_member(name: str, kind: str) -> pl.DataFrame:
    prefix, method = name.split("::", 1)
    table = pl.read_parquet(inp / cfg["sources"][prefix][kind])
    table = table.filter((pl.col("endpoint") == cfg["endpoint"]) & (pl.col("method") == method))
    column = prediction_column(table)
    if kind == "test":
        table = table.group_by("compound_id").agg(pl.col(column).mean().alias(name))
    else:
        table = table.select("compound_id", "fold", "y_true", pl.col(column).alias(name)).sort(
            "compound_id"
        )
    return table


members = sorted(weights)
oof = read_member(members[0], "oof")
test = read_member(members[0], "test")
for name in members[1:]:
    oof = oof.join(
        read_member(name, "oof").drop("fold", "y_true"),
        on="compound_id",
        validate="1:1",
    )
    test = test.join(read_member(name, "test"), on="compound_id", validate="1:1")

y = oof["y_true"].to_numpy().astype(int)
folds = oof["fold"].to_numpy().astype(int)
ensemble_oof = sum(weights[name] * oof[name].to_numpy() for name in members)
ensemble_test = sum(weights[name] * test[name].to_numpy() for name in members)
oof = oof.with_columns(pl.Series("caruana_auc", ensemble_oof))
test = test.with_columns(pl.Series("caruana_auc", ensemble_test))
all_names = members + ["caruana_auc"]


def fold_rank(values: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = rankdata(values[mask], method="average")
    return result


def calls(values: np.ndarray) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(positions))))
        order = np.argsort(-values[positions], kind="stable")
        result[positions[order[:count]]] = True
    return result


ranked = {name: fold_rank(oof[name].to_numpy()) for name in all_names}
called = {name: calls(oof[name].to_numpy()) for name in all_names}
metric_rows = []
for name in all_names:
    probability = oof[name].to_numpy()
    binary = called[name]
    metric_rows.append(
        {
            "member": name,
            "weight": weights.get(name, 1.0),
            "AUROC": float(roc_auc_score(y, probability)),
            "AP": float(average_precision_score(y, probability)),
            "MCC_at_fold_rate": float(matthews_corrcoef(y, binary)),
            "pearson_with_label": float(np.corrcoef(probability, y)[0, 1]),
            "spearman_with_ensemble": float(
                spearmanr(ranked[name], ranked["caruana_auc"]).statistic
            ),
        }
    )

pair_rows = []
for left_index, left in enumerate(all_names):
    for right in all_names[left_index + 1 :]:
        left_calls, right_calls = called[left], called[right]
        pair_rows.append(
            {
                "left": left,
                "right": right,
                "oof_pearson": float(np.corrcoef(oof[left], oof[right])[0, 1]),
                "oof_fold_rank_spearman": float(spearmanr(ranked[left], ranked[right]).statistic),
                "test_spearman": float(spearmanr(test[left], test[right]).statistic),
                "oof_call_jaccard": float(
                    np.logical_and(left_calls, right_calls).sum()
                    / np.logical_or(left_calls, right_calls).sum()
                ),
            }
        )

pl.DataFrame(metric_rows).sort("weight", descending=True).write_csv(out / "member-metrics.csv")
pl.DataFrame(pair_rows).write_csv(out / "pairwise-correlations.csv")
oof.select("compound_id", "fold", "y_true", *all_names).write_parquet(out / "member-oof.parquet")
test.select("compound_id", *all_names).write_parquet(out / "member-test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "endpoint": cfg["endpoint"],
            "recipe": cfg["recipe"],
            "weights": weights,
            "selection_steps": cfg["selection_steps"],
            "interpretation": (
                "Correlations describe saved predictions. Full-OOF deployment weights were "
                "selected on these OOF labels, so ensemble standalone metrics are descriptive, "
                "not an honest nested estimate."
            ),
        },
        indent=2,
    )
    + "\n"
)
print(pl.DataFrame(metric_rows).sort("weight", descending=True), flush=True)
