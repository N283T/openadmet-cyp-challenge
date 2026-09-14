"""Diverse non-LightGBM CYP2D6 rankers on the fixed UMAP folds."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from catboost import CatBoostClassifier
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import rankdata
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoint = "CYP2D6_is_TDI"


def read(name):
    return pl.read_parquet(inp / cfg[name])


def matrix(ids, table, prefix=None):
    columns = [
        column
        for column in table.columns
        if column != "compound_id" and (prefix is None or column.startswith(prefix))
    ]
    result = (
        ids.join(
            table.select("compound_id", *columns),
            on="compound_id",
            how="left",
            validate="1:1",
            maintain_order="left",
        )
        .select(columns)
        .to_numpy()
        .astype(np.float32)
    )
    result[~np.isfinite(result)] = np.nan
    return result


def unit_rank(values):
    return rankdata(values, method="average") / (len(values) + 1.0)


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


def fixed_rate_calls(values, folds, rate):
    calls = np.zeros(len(values), dtype=bool)
    for fold in range(5):
        positions = np.flatnonzero(folds == fold)
        count = int(round(rate * len(positions)))
        order = np.argsort(-values[positions], kind="stable")
        calls[positions[order[:count]]] = True
    return calls


def tanimoto_knn(query_fps, train_fps, train_y, neighbors):
    result = {k: np.empty(len(query_fps), dtype=float) for k in neighbors}
    for row, fingerprint in enumerate(query_fps):
        similarities = np.asarray(DataStructs.BulkTanimotoSimilarity(fingerprint, train_fps))
        order = np.argsort(-similarities, kind="stable")
        for k in neighbors:
            selected = order[:k]
            weights = similarities[selected] ** 2
            result[k][row] = float(np.average(train_y[selected], weights=weights + 1e-12))
    return result


compounds = read("compounds").sort("compound_id")
ids = compounds.select("compound_id")
labels = read("labels").filter(pl.col("endpoint") == endpoint)
folds = (
    read("assignments").filter(pl.col("strategy") == "umap").select("compound_id", "fold", "group")
)
data = labels.join(folds, on="compound_id", validate="1:1").sort("compound_id")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
index = {compound_id: row for row, compound_id in enumerate(ids["compound_id"])}
data_index = np.asarray([index[compound_id] for compound_id in data["compound_id"]])
test_index = np.asarray([index[compound_id] for compound_id in test["compound_id"]])
y = data["value"].to_numpy().astype(int)
fold_ids = data["fold"].to_numpy()
target_prior = json.loads((inp / cfg["lb_summary"]).read_text())["leaderboard_priors"][endpoint]

mordred = matrix(
    ids,
    pl.concat([pl.read_parquet(inp / filename) for filename in cfg["mordred"]]),
)
class_embedding = matrix(ids, read("class_emb"), "emb_")
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fingerprints = [
    generator.GetFingerprint(Chem.MolFromSmiles(s)) for s in compounds["standardized_smiles"]
]
morgan = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(s))
        for s in compounds["standardized_smiles"]
    ]
).astype(np.float32)
views = {
    "morgan": morgan,
    "mordred": mordred,
    "external_class": class_embedding,
    "fp_mordred": np.column_stack([morgan, mordred]),
}

oof_predictions = {}
test_predictions = {}
model_records = []
for fold in range(5):
    train = fold_ids != fold
    validation = ~train
    train_index = data_index[train]
    validation_index = data_index[validation]
    for view_name, features in views.items():
        model = CatBoostClassifier(
            iterations=350,
            depth=5,
            learning_rate=0.03,
            l2_leaf_reg=5.0,
            loss_function="Logloss",
            random_seed=20260913 + fold,
            thread_count=3,
            verbose=False,
            allow_writing_files=False,
        )
        model.fit(features[train_index], y[train])
        method = f"catboost_{view_name}"
        oof_predictions.setdefault(method, np.empty(len(y)))[validation] = model.predict_proba(
            features[validation_index]
        )[:, 1]
        test_predictions.setdefault(method, []).append(
            model.predict_proba(features[test_index])[:, 1]
        )
        model.save_model(out / f"{method}-fold{fold}.cbm")
        model_records.append({"method": method, "fold": fold, "dimensions": features.shape[1]})

    for leaf in [2, 5]:
        model = ExtraTreesClassifier(
            n_estimators=800,
            max_features="sqrt",
            min_samples_leaf=leaf,
            class_weight=None,
            n_jobs=3,
            random_state=20260913 + fold,
        ).fit(morgan[train_index], y[train])
        method = f"extra_morgan_leaf{leaf}"
        oof_predictions.setdefault(method, np.empty(len(y)))[validation] = model.predict_proba(
            morgan[validation_index]
        )[:, 1]
        test_predictions.setdefault(method, []).append(
            model.predict_proba(morgan[test_index])[:, 1]
        )

    for regularization in [0.05, 0.2, 1.0]:
        model = LogisticRegression(
            C=regularization,
            max_iter=3000,
            solver="liblinear",
            random_state=20260913 + fold,
        ).fit(morgan[train_index], y[train])
        method = f"logistic_morgan_c{regularization:g}"
        oof_predictions.setdefault(method, np.empty(len(y)))[validation] = model.predict_proba(
            morgan[validation_index]
        )[:, 1]
        test_predictions.setdefault(method, []).append(
            model.predict_proba(morgan[test_index])[:, 1]
        )

    train_fps = [fingerprints[i] for i in train_index]
    validation_scores = tanimoto_knn(
        [fingerprints[i] for i in validation_index], train_fps, y[train], [5, 15, 50]
    )
    test_scores = tanimoto_knn(
        [fingerprints[i] for i in test_index], train_fps, y[train], [5, 15, 50]
    )
    for neighbors in [5, 15, 50]:
        method = f"tanimoto_knn{neighbors}"
        oof_predictions.setdefault(method, np.empty(len(y)))[validation] = validation_scores[
            neighbors
        ]
        test_predictions.setdefault(method, []).append(test_scores[neighbors])
    print(f"fold {fold} complete", flush=True)

current_oof = (
    read("current_oof")
    .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == "all_plain_pfn_graph_mean"))
    .sort("compound_id")
)
assert current_oof["compound_id"].to_list() == data["compound_id"].to_list()
current_test_folds = read("current_test").filter(
    (pl.col("endpoint") == endpoint) & (pl.col("method") == "all_plain_pfn_graph_mean")
)
current_test = (
    current_test_folds.group_by("compound_id")
    .agg(pl.col("prediction").mean())
    .join(test, on="compound_id", validate="1:1", maintain_order="right")
)
assert current_test["compound_id"].to_list() == test["compound_id"].to_list()

# Rank within folds before blending so fold-specific probability scales do not
# leak into the comparison.
ranked_oof = {}
for method, prediction in {
    "current": current_oof["prediction"].to_numpy(),
    **oof_predictions,
}.items():
    ranked = np.empty(len(y))
    for fold in range(5):
        mask = fold_ids == fold
        ranked[mask] = unit_rank(prediction[mask])
    ranked_oof[method] = ranked
ranked_test = {"current": unit_rank(current_test["prediction"].to_numpy())}
for method, predictions in test_predictions.items():
    ranked_test[method] = unit_rank(np.mean(predictions, axis=0))

for method in list(oof_predictions):
    ranked_oof[f"current_{method}_rank_mean"] = (ranked_oof["current"] + ranked_oof[method]) / 2
    ranked_test[f"current_{method}_rank_mean"] = (ranked_test["current"] + ranked_test[method]) / 2
    ranked_oof[f"current2_{method}_rank_mean"] = (
        2 * ranked_oof["current"] + ranked_oof[method]
    ) / 3
    ranked_test[f"current2_{method}_rank_mean"] = (
        2 * ranked_test["current"] + ranked_test[method]
    ) / 3

rates = {
    "current_full_rate": 43 / 750,
    "top10_median_rate": cfg["top10_median_rate"],
    "target_prior": target_prior,
}
metric_rows = []
for method, prediction in ranked_oof.items():
    for rate_name, rate in rates.items():
        binary = fixed_rate_calls(prediction, fold_ids, rate)
        metric_rows.append(
            {
                "method": method,
                "rate_policy": rate_name,
                "rate": rate,
                "source_MCC": float(matthews_corrcoef(y, binary)),
                "expected_target_MCC": expected_mcc(y, binary, target_prior),
                "AUROC": float(roc_auc_score(y, prediction)),
                "AP": float(average_precision_score(y, prediction)),
                "oof_positives": int(binary.sum()),
            }
        )
metrics = pl.DataFrame(metric_rows).sort(
    "rate_policy", "expected_target_MCC", descending=[False, True]
)
metrics.write_csv(out / "metrics.csv")

oof_rows = []
test_rows = []
for method, prediction in ranked_oof.items():
    oof_rows.append(
        data.select("compound_id", "fold", "group", pl.col("value").alias("y_true")).with_columns(
            pl.lit(endpoint).alias("endpoint"),
            pl.lit(method).alias("method"),
            pl.Series("prediction", prediction),
        )
    )
for method, prediction in ranked_test.items():
    test_rows.append(
        test.with_columns(
            pl.lit(endpoint).alias("endpoint"),
            pl.lit(method).alias("method"),
            pl.Series("prediction", prediction),
        )
    )
pl.concat(oof_rows).write_parquet(out / "oof.parquet")
pl.concat(test_rows).write_parquet(out / "test.parquet")
(out / "models.json").write_text(json.dumps(model_records, indent=2) + "\n")
(out / "summary.json").write_text(
    json.dumps(
        {
            "endpoint": endpoint,
            "split_run": cfg["split_run"],
            "target_prior": target_prior,
            "base_methods": sorted(oof_predictions),
            "rates": rates,
            "selection_warning": (
                "Model/rate screen on the same five outer-fold predictions; "
                "projected target MCC assumes transferable TPR/FPR."
            ),
        },
        indent=2,
    )
    + "\n"
)
print(metrics.group_by("rate_policy").head(8).write_csv(), flush=True)
