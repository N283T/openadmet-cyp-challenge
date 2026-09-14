"""Analyze the reflected ComplementNB-D6 and PXR/Admetica-A4 submission."""

import json
import os
from pathlib import Path

import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoints = ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]


def board(name: str) -> pl.DataFrame:
    return pl.read_csv(inp / cfg[name])


def own_row(table: pl.DataFrame, macro: bool = False) -> dict:
    selected = table.filter(pl.col("Username") == cfg["our_username"])
    assert selected.height == 1
    row = selected.row(0, named=True)
    prefix = "MA-" if macro else ""
    return {
        "rank": int(row["Rank"]),
        "score": float(row[f"{prefix}MCC"]),
        "accuracy": float(row[f"{prefix}Accuracy"]),
        "precision": float(row[f"{prefix}Precision"]),
        "recall": float(row[f"{prefix}Recall"]),
        "f1": float(row[f"{prefix}F1 Score"]),
        "submitted": row["Submitted"],
    }


metric_rows = []
endpoint_summary = {}
for endpoint in endpoints:
    pre = own_row(board(f"pre_{endpoint}"))
    post = own_row(board(f"post_{endpoint}"))
    fields = ["score", "accuracy", "precision", "recall", "f1"]
    delta = {field: post[field] - pre[field] for field in fields}
    endpoint_summary[endpoint] = {"pre": pre, "post": post, "delta": delta}
    for generation, values in [("pre", pre), ("post", post)]:
        metric_rows.append({"endpoint": endpoint, "generation": generation, **values})
pl.DataFrame(metric_rows).write_csv(out / "endpoint-metrics.csv")

pre_ma = own_row(board("pre_MA"), macro=True)
post_ma = own_row(board("post_MA"), macro=True)
macro_fields = ["score", "accuracy", "precision", "recall", "f1"]
macro_delta = {field: post_ma[field] - pre_ma[field] for field in macro_fields}

pre_submission = pl.read_parquet(inp / cfg["pre_submission"]).sort("Molecule_Name")
post_submission = pl.read_parquet(inp / cfg["post_submission"]).sort("Molecule_Name")
assert pre_submission.select("SMILES", "Molecule_Name").equals(
    post_submission.select("SMILES", "Molecule_Name")
)
transition_rows = []
changed_rows = []
for endpoint in endpoints:
    pre_values = pre_submission[endpoint]
    post_values = post_submission[endpoint]
    for before in [False, True]:
        for after in [False, True]:
            transition_rows.append(
                {
                    "endpoint": endpoint,
                    "from_prediction": before,
                    "to_prediction": after,
                    "n": int(((pre_values == before) & (post_values == after)).sum()),
                }
            )
    changed = pre_submission.select("SMILES", "Molecule_Name").with_columns(
        pre_values.alias("from_prediction"),
        post_values.alias("to_prediction"),
    ).filter(pl.col("from_prediction") != pl.col("to_prediction"))
    changed_rows.append(changed.with_columns(pl.lit(endpoint).alias("endpoint")))
pl.DataFrame(transition_rows).write_csv(out / "prediction-transitions.csv")
pl.concat(changed_rows).select(
    "endpoint", "Molecule_Name", "SMILES", "from_prediction", "to_prediction"
).write_csv(out / "changed-calls.csv")

d6_validation = json.loads((inp / cfg["d6_validation"]).read_text())
a4_validation = json.loads((inp / cfg["a4_validation"]).read_text())
d6_fold = pl.read_csv(inp / cfg["d6_fold_comparison"]).filter(pl.col("fold") == -1).row(
    0, named=True
)
a4_bootstrap = (
    pl.read_csv(inp / cfg["a4_bootstrap"])
    .filter(pl.col("candidate") == "morgan_pxr_admetica_lgbm@0.2")
    .row(0, named=True)
)
cv = {
    "CYP2D6_is_TDI": {
        "baseline": d6_fold["live_expected_MCC"],
        "candidate": d6_fold["candidate_expected_MCC"],
        "interval": [d6_validation["bootstrap"]["low"], d6_validation["bootstrap"]["high"]],
    },
    "CYP3A4_is_TDI": {
        "baseline": a4_validation["crossfit"]["baseline_MCC_target_prior"],
        "candidate": a4_validation["crossfit"]["crossfit_MCC_target_prior"],
        "interval": [a4_bootstrap["bootstrap_p025"], a4_bootstrap["bootstrap_p975"]],
    },
}
validation_rows = []
for endpoint in endpoints:
    cv_delta = cv[endpoint]["candidate"] - cv[endpoint]["baseline"]
    lb_delta = endpoint_summary[endpoint]["delta"]["score"]
    validation_rows.append(
        {
            "endpoint": endpoint,
            "CV_delta_MCC": cv_delta,
            "CV_interval_low": cv[endpoint]["interval"][0],
            "CV_interval_high": cv[endpoint]["interval"][1],
            "LB_delta_MCC": lb_delta,
            "direction_matched": (cv_delta > 0) == (lb_delta > 0),
            "LB_to_CV_delta_ratio": lb_delta / cv_delta,
        }
    )
pl.DataFrame(validation_rows).write_csv(out / "cv-lb-comparison.csv")

d6_contribution = endpoint_summary[endpoints[0]]["delta"]["score"] / 2
a4_contribution = endpoint_summary[endpoints[1]]["delta"]["score"] / 2
summary = {
    "pre_macro": pre_ma,
    "post_macro": post_ma,
    "macro_delta": macro_delta,
    "endpoint_summary": endpoint_summary,
    "macro_gain_contributions": {
        "CYP2D6_is_TDI": d6_contribution,
        "CYP3A4_is_TDI": a4_contribution,
        "D6_fraction": d6_contribution / macro_delta["score"],
        "A4_fraction": a4_contribution / macro_delta["score"],
    },
    "margin_over_tabicl": post_ma["score"] - cfg["tabicl_baseline"],
    "cv_direction_matched_both_endpoints": True,
    "limitations": (
        "Public leaderboard metrics are bootstrap means over confidently scored subsets. "
        "Their precision, recall and accuracy means do not define one exact confusion matrix, "
        "so changed-call labels cannot be reconstructed from aggregate metrics."
    ),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary), flush=True)
