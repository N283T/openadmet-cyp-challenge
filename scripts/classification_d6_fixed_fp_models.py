"""Compare classical model families on one fixed Morgan CYP2D6 representation."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import rankdata
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression, RidgeClassifier, SGDClassifier
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from sklearn.naive_bayes import BernoulliNB, ComplementNB
from sklearn.svm import SVC, LinearSVC

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoint = "CYP2D6_is_TDI"
rate = 49 / 750


def read(name):
    return pl.read_parquet(inp / cfg[name])


def unit_rank(values):
    return rankdata(values, method="average") / (len(values) + 1.0)


def score(model, features):
    if hasattr(model, "predict_proba"):
        return model.predict_proba(features)[:, 1]
    return model.decision_function(features)


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


def fixed_rate_calls(values, folds):
    calls = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = int(round(rate * len(positions)))
        order = np.argsort(-values[positions], kind="stable")
        calls[positions[order[:count]]] = True
    return calls


def tanimoto_kernel(left, right):
    intersection = left @ right.T
    denominator = left.sum(axis=1)[:, None] + right.sum(axis=1)[None, :] - intersection
    return np.divide(
        intersection,
        denominator,
        out=np.zeros_like(intersection, dtype=np.float32),
        where=denominator > 0,
    )


compounds = read("compounds").sort("compound_id")
labels = read("labels").filter(pl.col("endpoint") == endpoint)
fold_table = (
    read("assignments").filter(pl.col("strategy") == "umap").select("compound_id", "fold", "group")
)
data = labels.join(fold_table, on="compound_id", validate="1:1").sort("compound_id")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
index = {compound_id: row for row, compound_id in enumerate(compounds["compound_id"])}
data_index = np.asarray([index[compound_id] for compound_id in data["compound_id"]])
test_index = np.asarray([index[compound_id] for compound_id in test["compound_id"]])
y = data["value"].to_numpy().astype(int)
folds = data["fold"].to_numpy()
target_prior = cfg["target_prior"]

generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
all_features = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in compounds["standardized_smiles"]
    ]
).astype(np.float32)
x = all_features[data_index]
xtest = all_features[test_index]
combined = np.vstack([x, xtest])
tanimoto = tanimoto_kernel(combined, combined)

oof_predictions = {}
test_predictions = {}
model_metadata = []
for fold in range(5):
    train = folds != fold
    validation = ~train
    models = []
    for c in [0.01, 0.05, 0.2, 1.0]:
        models.append(
            (
                f"logistic_c{c:g}",
                LogisticRegression(
                    C=c,
                    max_iter=4000,
                    solver="liblinear",
                    random_state=20260913 + fold,
                ),
            )
        )
    for c in [0.001, 0.01, 0.1, 1.0]:
        models.append((f"linear_svm_c{c:g}", LinearSVC(C=c, dual="auto")))
    for c in [0.01, 0.1]:
        models.append(
            (
                f"linear_svm_balanced_c{c:g}",
                LinearSVC(C=c, dual="auto", class_weight="balanced"),
            )
        )
    for c in [0.1, 1.0, 10.0]:
        models.append((f"rbf_svm_c{c:g}", SVC(C=c, kernel="rbf")))
    models.append(("rbf_svm_balanced_c1", SVC(C=1.0, kernel="rbf", class_weight="balanced")))
    for leaf in [1, 3, 8]:
        for balanced in [False, True]:
            suffix = "_balanced" if balanced else ""
            common = {
                "n_estimators": 500,
                "max_features": "sqrt",
                "min_samples_leaf": leaf,
                "class_weight": "balanced" if balanced else None,
                "n_jobs": 3,
                "random_state": 20260913 + fold,
            }
            models.append((f"random_forest_leaf{leaf}{suffix}", RandomForestClassifier(**common)))
            models.append((f"extra_trees_leaf{leaf}{suffix}", ExtraTreesClassifier(**common)))
    for alpha in [0.1, 1.0, 10.0]:
        models.append((f"bernoulli_nb_a{alpha:g}", BernoulliNB(alpha=alpha)))
        models.append((f"complement_nb_a{alpha:g}", ComplementNB(alpha=alpha)))
    for alpha in [1.0, 10.0, 100.0]:
        models.append((f"ridge_a{alpha:g}", RidgeClassifier(alpha=alpha)))
    models.append(("ridge_balanced_a10", RidgeClassifier(alpha=10, class_weight="balanced")))
    for loss in ["log_loss", "modified_huber"]:
        for alpha in [1e-5, 1e-4, 1e-3]:
            models.append(
                (
                    f"sgd_{loss}_a{alpha:g}",
                    SGDClassifier(
                        loss=loss,
                        alpha=alpha,
                        max_iter=5000,
                        tol=1e-5,
                        random_state=20260913 + fold,
                    ),
                )
            )
    models.append(
        (
            "lda_shrinkage_auto",
            LinearDiscriminantAnalysis(solver="lsqr", shrinkage="auto"),
        )
    )
    for name, model in models:
        model.fit(x[train], y[train])
        oof_predictions.setdefault(name, np.empty(len(y)))[validation] = score(model, x[validation])
        test_predictions.setdefault(name, []).append(score(model, xtest))
        joblib.dump(model, out / f"{name}-fold{fold}.joblib", compress=3)
        model_metadata.append({"method": name, "fold": fold, "kind": type(model).__name__})

    train_positions = np.flatnonzero(train)
    validation_positions = np.flatnonzero(validation)
    test_positions = np.arange(len(y), len(y) + len(test))
    for c in [0.1, 1.0, 10.0]:
        for balanced in [False, True] if c == 1.0 else [False]:
            suffix = "_balanced" if balanced else ""
            name = f"tanimoto_svm_c{c:g}{suffix}"
            model = SVC(
                C=c,
                kernel="precomputed",
                class_weight="balanced" if balanced else None,
            ).fit(tanimoto[np.ix_(train_positions, train_positions)], y[train])
            oof_predictions.setdefault(name, np.empty(len(y)))[validation] = (
                model.decision_function(tanimoto[np.ix_(validation_positions, train_positions)])
            )
            test_predictions.setdefault(name, []).append(
                model.decision_function(tanimoto[np.ix_(test_positions, train_positions)])
            )
            joblib.dump(model, out / f"{name}-fold{fold}.joblib", compress=3)
            model_metadata.append({"method": name, "fold": fold, "kind": "SVC-tanimoto"})
    print(f"fold {fold} complete", flush=True)

original_oof = (
    read("original_oof")
    .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == "all_plain_pfn_graph_mean"))
    .sort("compound_id")
)
original_test = (
    read("original_test")
    .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == "all_plain_pfn_graph_mean"))
    .group_by("compound_id")
    .agg(pl.col("prediction").mean())
    .join(test, on="compound_id", validate="1:1", maintain_order="right")
)
live_method = cfg["live_method"]
live_oof = read("live_oof").filter(pl.col("method") == live_method).sort("compound_id")
live_test = read("live_test").filter(pl.col("method") == live_method).sort("compound_id")
assert original_oof["compound_id"].to_list() == data["compound_id"].to_list()
assert live_oof["compound_id"].to_list() == data["compound_id"].to_list()

ranked_oof = {
    "original": np.empty(len(y)),
    "live": live_oof["prediction"].to_numpy(),
}
for fold in range(5):
    mask = folds == fold
    ranked_oof["original"][mask] = unit_rank(
        original_oof.filter(pl.col("fold") == fold)["prediction"]
    )
ranked_test = {
    "original": unit_rank(original_test["prediction"].to_numpy()),
    "live": live_test["prediction"].to_numpy(),
}
for method, prediction in oof_predictions.items():
    ranked = np.empty(len(y))
    for fold in range(5):
        mask = folds == fold
        ranked[mask] = unit_rank(prediction[mask])
    ranked_oof[method] = ranked
    fold_ranks = np.stack([unit_rank(values) for values in test_predictions[method]])
    ranked_test[method] = unit_rank(fold_ranks.mean(axis=0))
    ranked_oof[f"original2_{method}"] = (2 * ranked_oof["original"] + ranked) / 3
    ranked_test[f"original2_{method}"] = (2 * ranked_test["original"] + ranked_test[method]) / 3
    ranked_oof[f"live2_{method}"] = (2 * ranked_oof["live"] + ranked) / 3
    ranked_test[f"live2_{method}"] = (2 * ranked_test["live"] + ranked_test[method]) / 3

screen_rows = []
for method, prediction in ranked_oof.items():
    calls = fixed_rate_calls(prediction, folds)
    screen_rows.append(
        {
            "method": method,
            "source_MCC": float(matthews_corrcoef(y, calls)),
            "expected_target_MCC": expected_mcc(y, calls, target_prior),
            "AUROC": float(roc_auc_score(y, prediction)),
            "AP": float(average_precision_score(y, prediction)),
        }
    )
screen = pl.DataFrame(screen_rows).sort("expected_target_MCC", descending=True)
screen.write_csv(out / "screen.csv")

incremental = ["live", *[method for method in ranked_oof if method.startswith("live2_")]]
selection_rows = []
deployment_rows = []
for top_k in [1, 3, 5, 10]:
    crossfit_prediction = np.empty(len(y))
    fold_choices = []
    for heldout in range(5):
        selection_mask = folds != heldout
        scored = []
        for method in incremental:
            calls = fixed_rate_calls(ranked_oof[method][selection_mask], folds[selection_mask])
            scored.append((expected_mcc(y[selection_mask], calls, target_prior), method))
        chosen = [method for _, method in sorted(scored, reverse=True)[:top_k]]
        fold_choices.append(chosen)
        mask = folds == heldout
        crossfit_prediction[mask] = np.mean([ranked_oof[method][mask] for method in chosen], axis=0)
    calls = fixed_rate_calls(crossfit_prediction, folds)
    full_scored = []
    for method in incremental:
        method_calls = fixed_rate_calls(ranked_oof[method], folds)
        full_scored.append((expected_mcc(y, method_calls, target_prior), method))
    deployment_methods = [method for _, method in sorted(full_scored, reverse=True)[:top_k]]
    test_prediction = np.mean([ranked_test[method] for method in deployment_methods], axis=0)
    test_calls = np.zeros(len(test), dtype=bool)
    test_calls[np.argsort(-test_prediction, kind="stable")[:49]] = True
    variant = f"top{top_k}"
    selection_rows.append(
        {
            "variant": variant,
            "crossfit_source_MCC": float(matthews_corrcoef(y, calls)),
            "crossfit_expected_target_MCC": expected_mcc(y, calls, target_prior),
            "crossfit_AUROC": float(roc_auc_score(y, crossfit_prediction)),
            "deployment_methods": json.dumps(deployment_methods),
            "fold_methods": json.dumps(fold_choices),
        }
    )
    deployment_rows.append(
        test.with_columns(
            pl.lit(endpoint).alias("endpoint"),
            pl.lit(variant).alias("method"),
            pl.Series("prediction", test_prediction),
            pl.Series("binary", test_calls),
        )
    )
selection = pl.DataFrame(selection_rows).sort("crossfit_expected_target_MCC", descending=True)
selection.write_csv(out / "crossfit-selection.csv")
pl.concat(deployment_rows).write_parquet(out / "deployment.parquet")

oof_rows = []
for method, prediction in ranked_oof.items():
    oof_rows.append(
        data.select("compound_id", "fold", "group", pl.col("value").alias("y_true")).with_columns(
            pl.lit(endpoint).alias("endpoint"),
            pl.lit(method).alias("method"),
            pl.Series("prediction", prediction),
        )
    )
pl.concat(oof_rows).write_parquet(out / "oof.parquet")
(out / "models.json").write_text(json.dumps(model_metadata, indent=2) + "\n")
(out / "summary.json").write_text(
    json.dumps(
        {
            "feature": "Morgan radius2 2048-bit",
            "base_models": len(oof_predictions),
            "evaluated_ranks": len(ranked_oof),
            "positive_calls": 49,
            "selection": "Each heldout fold selects only from the other four OOF folds.",
        },
        indent=2,
    )
    + "\n"
)
print(selection.write_csv(), flush=True)
