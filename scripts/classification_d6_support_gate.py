"""Evaluate support-aware routing among saved CYP2D6 predictions."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoint = "CYP2D6_is_TDI"
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def read(name):
    return pl.read_parquet(inp / cfg[name])


def unit_rank(values):
    return rankdata(values, method="average") / (len(values) + 1.0)


def fold_rank(values, folds):
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = unit_rank(values[mask])
    return result


def fixed_rate_calls(values, folds, rate):
    calls = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = int(round(rate * len(positions)))
        order = np.argsort(-values[positions], kind="stable")
        calls[positions[order[:count]]] = True
    return calls


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


def table_prediction(table, method, ids):
    selected = table.filter(pl.col("method") == method).sort("compound_id")
    assert selected["compound_id"].to_list() == ids["compound_id"].to_list(), method
    return selected["prediction"].to_numpy()


def model_score(model, features):
    if hasattr(model, "predict_proba"):
        return model.predict_proba(features)[:, 1]
    return model.decision_function(features)


compounds = read("compounds").sort("compound_id")
all_ids = compounds.select("compound_id")
labels = read("labels").filter(pl.col("endpoint") == endpoint)
assignments = read("assignments").filter(pl.col("strategy") == "umap")
data = labels.join(
    assignments.select("compound_id", "fold", "group"), on="compound_id", validate="1:1"
).sort("compound_id")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
id_to_index = {
    compound_id: row for row, compound_id in enumerate(all_ids["compound_id"])
}
data_index = np.asarray([id_to_index[value] for value in data["compound_id"]])
test_index = np.asarray([id_to_index[value] for value in test["compound_id"]])
y = data["value"].to_numpy().astype(int)
folds = data["fold"].to_numpy().astype(int)

molecules = [Chem.MolFromSmiles(value) for value in compounds["standardized_smiles"]]
fingerprints = [generator.GetFingerprint(molecule) for molecule in molecules]
x = np.stack([generator.GetFingerprintAsNumPy(molecule) for molecule in molecules]).astype(
    np.float32
)

support = np.empty(len(data), dtype=float)
for fold in range(5):
    train = folds != fold
    validation = ~train
    train_fps = [fingerprints[data_index[row]] for row in np.flatnonzero(train)]
    for row in np.flatnonzero(validation):
        support[row] = max(
            DataStructs.BulkTanimotoSimilarity(fingerprints[data_index[row]], train_fps)
        )
support_rank = fold_rank(support, folds)

full_train_fps = [fingerprints[index] for index in data_index]
test_support = np.asarray(
    [
        max(DataStructs.BulkTanimotoSimilarity(fingerprints[index], full_train_fps))
        for index in test_index
    ]
)
test_support_rank = unit_rank(test_support)

diverse_oof = read("diverse_oof")
diverse_test = read("diverse_test")
fixed_oof = read("fixed_oof")
fp_oof = read("fp_oof")
fp_test = read("fp_test")

oof_predictions = {
    "live": table_prediction(diverse_oof, cfg["methods"]["live"], data.select("compound_id")),
    "original": table_prediction(
        diverse_oof, cfg["methods"]["original"], data.select("compound_id")
    ),
    "tanimoto_knn5": table_prediction(
        diverse_oof, cfg["methods"]["tanimoto_knn5"], data.select("compound_id")
    ),
    "catboost_mordred": table_prediction(
        diverse_oof, cfg["methods"]["catboost_mordred"], data.select("compound_id")
    ),
    "catboost_external": table_prediction(
        diverse_oof, cfg["methods"]["catboost_external"], data.select("compound_id")
    ),
    "complement_blend": table_prediction(
        fixed_oof, cfg["methods"]["complement_blend"], data.select("compound_id")
    ),
    "rbf_blend": table_prediction(
        fixed_oof, cfg["methods"]["rbf_blend"], data.select("compound_id")
    ),
    "r3_count": table_prediction(
        fp_oof, cfg["methods"]["r3_count"], data.select("compound_id")
    ),
    "chemprop_mixed": table_prediction(
        fp_oof, cfg["methods"]["chemprop_mixed"], data.select("compound_id")
    ),
    "r3_chemprop_fusion": table_prediction(
        fp_oof, cfg["methods"]["r3_chemprop_fusion"], data.select("compound_id")
    ),
}

test_predictions = {
    "live": table_prediction(diverse_test, cfg["methods"]["live"], test),
    "original": table_prediction(diverse_test, cfg["methods"]["original"], test),
    "tanimoto_knn5": table_prediction(diverse_test, cfg["methods"]["tanimoto_knn5"], test),
    "catboost_mordred": table_prediction(
        diverse_test, cfg["methods"]["catboost_mordred"], test
    ),
    "catboost_external": table_prediction(
        diverse_test, cfg["methods"]["catboost_external"], test
    ),
    "r3_count": table_prediction(fp_test, cfg["methods"]["r3_count"], test),
    "chemprop_mixed": table_prediction(fp_test, cfg["methods"]["chemprop_mixed"], test),
    "r3_chemprop_fusion": table_prediction(
        fp_test, cfg["methods"]["r3_chemprop_fusion"], test
    ),
}

# Reconstruct the two fixed-feature test blends, whose source run retained models
# and deployment predictions but did not write every raw test rank.
base_test = {}
for family in ["complement", "rbf"]:
    fold_predictions = []
    for filename in cfg["fixed_models"][family]:
        model = joblib.load(inp / filename)
        fold_predictions.append(unit_rank(model_score(model, x[test_index])))
    base_test[family] = unit_rank(np.mean(fold_predictions, axis=0))
test_predictions["complement_blend"] = (
    2 * test_predictions["live"] + base_test["complement"]
) / 3
test_predictions["rbf_blend"] = (2 * test_predictions["original"] + base_test["rbf"]) / 3

for name in oof_predictions:
    oof_predictions[name] = fold_rank(oof_predictions[name], folds)
    test_predictions[name] = unit_rank(test_predictions[name])

rate = cfg["positive_calls"] / len(test)
stratum_rows = []
bin_edges = cfg["support_bins"]
for name, prediction in oof_predictions.items():
    calls = fixed_rate_calls(prediction, folds, rate)
    for low, high in zip(bin_edges[:-1], bin_edges[1:], strict=True):
        mask = (support >= low) & (support < high)
        if mask.sum() == 0:
            continue
        row = {
            "method": name,
            "support_low": low,
            "support_high": high,
            "n": int(mask.sum()),
            "positives": int(y[mask].sum()),
            "positive_rate": float(y[mask].mean()),
            "global_calls_in_stratum": int(calls[mask].sum()),
            "global_call_precision_in_stratum": float(y[mask & calls].mean())
            if (mask & calls).sum()
            else None,
            "AUROC": float(roc_auc_score(y[mask], prediction[mask]))
            if len(np.unique(y[mask])) == 2
            else None,
            "AP": float(average_precision_score(y[mask], prediction[mask]))
            if y[mask].sum()
            else None,
        }
        stratum_rows.append(row)
pl.DataFrame(stratum_rows).write_csv(out / "support-stratum-metrics.csv")

low_methods = cfg["low_support_methods"]
high_methods = cfg["high_support_methods"]
gate_specs = [
    (quantile, low_method, high_method)
    for quantile in cfg["support_quantiles"]
    for low_method in low_methods
    for high_method in high_methods
]

selection_rows = []
deployment_rows = []
for top_k in cfg["top_k"]:
    crossfit_prediction = np.empty(len(y), dtype=float)
    fold_choices = []
    for heldout in range(5):
        selection = folds != heldout
        candidates = []
        for quantile, low_method, high_method in gate_specs:
            gated = np.where(
                support_rank < quantile,
                oof_predictions[low_method],
                oof_predictions[high_method],
            )
            calls = fixed_rate_calls(gated[selection], folds[selection], rate)
            score = expected_mcc(y[selection], calls, cfg["target_prior"])
            candidates.append((score, quantile, low_method, high_method))
        chosen = sorted(candidates, reverse=True)[:top_k]
        fold_choices.append(chosen)
        heldout_mask = folds == heldout
        fold_predictions = []
        for _, quantile, low_method, high_method in chosen:
            fold_predictions.append(
                np.where(
                    support_rank[heldout_mask] < quantile,
                    oof_predictions[low_method][heldout_mask],
                    oof_predictions[high_method][heldout_mask],
                )
            )
        crossfit_prediction[heldout_mask] = np.mean(fold_predictions, axis=0)

    crossfit_calls = fixed_rate_calls(crossfit_prediction, folds, rate)
    full_candidates = []
    for quantile, low_method, high_method in gate_specs:
        gated = np.where(
            support_rank < quantile,
            oof_predictions[low_method],
            oof_predictions[high_method],
        )
        calls = fixed_rate_calls(gated, folds, rate)
        full_candidates.append(
            (
                expected_mcc(y, calls, cfg["target_prior"]),
                quantile,
                low_method,
                high_method,
            )
        )
    deployment = sorted(full_candidates, reverse=True)[:top_k]
    deployment_predictions = []
    for _, quantile, low_method, high_method in deployment:
        deployment_predictions.append(
            np.where(
                test_support_rank < quantile,
                test_predictions[low_method],
                test_predictions[high_method],
            )
        )
    test_prediction = np.mean(deployment_predictions, axis=0)
    test_calls = np.zeros(len(test), dtype=bool)
    test_calls[np.argsort(-test_prediction, kind="stable")[: cfg["positive_calls"]]] = True
    variant = f"support_gate_top{top_k}"
    selection_rows.append(
        {
            "variant": variant,
            "crossfit_source_MCC": float(matthews_corrcoef(y, crossfit_calls)),
            "crossfit_expected_target_MCC": expected_mcc(
                y, crossfit_calls, cfg["target_prior"]
            ),
            "crossfit_AUROC": float(roc_auc_score(y, crossfit_prediction)),
            "crossfit_AP": float(average_precision_score(y, crossfit_prediction)),
            "fold_choices": json.dumps(fold_choices),
            "deployment": json.dumps(deployment),
        }
    )
    deployment_rows.append(
        test.with_columns(
            pl.lit(endpoint).alias("endpoint"),
            pl.lit(variant).alias("method"),
            pl.Series("prediction", test_prediction),
            pl.Series("binary", test_calls),
            pl.Series("nearest_train_tanimoto", test_support),
            pl.Series("support_rank", test_support_rank),
        )
    )

live_calls = fixed_rate_calls(oof_predictions["live"], folds, rate)
selection_rows.append(
    {
        "variant": "live_reference",
        "crossfit_source_MCC": float(matthews_corrcoef(y, live_calls)),
        "crossfit_expected_target_MCC": expected_mcc(y, live_calls, cfg["target_prior"]),
        "crossfit_AUROC": float(roc_auc_score(y, oof_predictions["live"])),
        "crossfit_AP": float(average_precision_score(y, oof_predictions["live"])),
        "fold_choices": None,
        "deployment": None,
    }
)
selection = pl.DataFrame(selection_rows).sort(
    "crossfit_expected_target_MCC", descending=True
)
selection.write_csv(out / "gate-selection.csv")
pl.concat(deployment_rows).write_parquet(out / "deployment.parquet")

test_summary = pl.DataFrame(
    {
        "set": ["oof", "test"],
        "n": [len(support), len(test_support)],
        "support_mean": [float(support.mean()), float(test_support.mean())],
        "support_median": [float(np.median(support)), float(np.median(test_support))],
        "below_0_3": [float((support < 0.3).mean()), float((test_support < 0.3).mean())],
        "below_0_4": [float((support < 0.4).mean()), float((test_support < 0.4).mean())],
    }
)
test_summary.write_csv(out / "support-distribution.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "endpoint": endpoint,
            "methods": list(oof_predictions),
            "gate_candidates": len(gate_specs),
            "selection": "Each heldout fold selects gates using the other four OOF folds.",
            "support": (
                "Maximum Morgan Tanimoto to fold-training molecules; gating uses "
                "within-fold percentile ranks to transfer to full-train test support."
            ),
            "uploaded": False,
        },
        indent=2,
    )
    + "\n"
)
print(test_summary.write_csv(), flush=True)
print(selection.write_csv(), flush=True)
