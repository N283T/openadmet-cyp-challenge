"""Classical CYP2D6 screen across fingerprints and frozen Chemprop embeddings."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from sklearn.preprocessing import StandardScaler

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
    values = (
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
    values[~np.isfinite(values)] = np.nan
    return values


def fingerprint_matrix(molecules, generator, count=False):
    function = generator.GetCountFingerprintAsNumPy if count else generator.GetFingerprintAsNumPy
    return np.stack([function(molecule) for molecule in molecules]).astype(np.float32)


def maccs_matrix(molecules):
    values = []
    for molecule in molecules:
        fingerprint = MACCSkeys.GenMACCSKeys(molecule)
        row = np.zeros(fingerprint.GetNumBits(), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fingerprint, row)
        values.append(row)
    return np.stack(values).astype(np.float32)


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
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = int(round(rate * len(positions)))
        order = np.argsort(-values[positions], kind="stable")
        calls[positions[order[:count]]] = True
    return calls


def fit_logistic(name, x_train, y_train, x_val, x_test, fold, c=0.05):
    model = LogisticRegression(
        C=c,
        max_iter=4000,
        solver="liblinear",
        random_state=20260913 + fold,
    ).fit(x_train, y_train)
    np.savez_compressed(
        out / f"{name}-fold{fold}.npz",
        coef=model.coef_,
        intercept=model.intercept_,
        classes=model.classes_,
    )
    return model.predict_proba(x_val)[:, 1], model.predict_proba(x_test)[:, 1]


compounds = read("compounds").sort("compound_id")
ids = compounds.select("compound_id")
labels = read("labels").filter(pl.col("endpoint") == endpoint)
fold_table = (
    read("assignments").filter(pl.col("strategy") == "umap").select("compound_id", "fold", "group")
)
data = labels.join(fold_table, on="compound_id", validate="1:1").sort("compound_id")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
index = {compound_id: row for row, compound_id in enumerate(ids["compound_id"])}
data_index = np.asarray([index[compound_id] for compound_id in data["compound_id"]])
test_index = np.asarray([index[compound_id] for compound_id in test["compound_id"]])
y = data["value"].to_numpy().astype(int)
folds = data["fold"].to_numpy()
target_prior = cfg["target_prior"]

molecules = [Chem.MolFromSmiles(smiles) for smiles in compounds["standardized_smiles"]]
fingerprints = {
    "morgan_r1_bit": fingerprint_matrix(
        molecules, rdFingerprintGenerator.GetMorganGenerator(radius=1, fpSize=2048)
    ),
    "morgan_r2_bit": fingerprint_matrix(
        molecules, rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    ),
    "morgan_r3_bit": fingerprint_matrix(
        molecules, rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
    ),
    "morgan_r2_chiral": fingerprint_matrix(
        molecules,
        rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True),
    ),
    "morgan_r2_count": fingerprint_matrix(
        molecules,
        rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048),
        count=True,
    ),
    "morgan_r3_count": fingerprint_matrix(
        molecules,
        rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048),
        count=True,
    ),
    "rdkit_path5": fingerprint_matrix(
        molecules, rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5, fpSize=2048)
    ),
    "rdkit_path7": fingerprint_matrix(
        molecules, rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=7, fpSize=2048)
    ),
    "atom_pair_count": fingerprint_matrix(
        molecules,
        rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048),
        count=True,
    ),
    "torsion_count": fingerprint_matrix(
        molecules,
        rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048),
        count=True,
    ),
    "maccs": maccs_matrix(molecules),
}
embeddings = {
    "chemprop_class": matrix(ids, read("class_emb"), "emb_"),
    "chemprop_mixed": matrix(ids, read("mixed_emb"), "emb_"),
}

oof_predictions = {}
test_predictions = {}
for fold in range(5):
    train = folds != fold
    validation = ~train
    train_index = data_index[train]
    validation_index = data_index[validation]
    for fp_name, features in fingerprints.items():
        for c in [0.02, 0.05, 0.2]:
            method = f"logistic_{fp_name}_c{c:g}"
            val_score, test_score = fit_logistic(
                method,
                features[train_index],
                y[train],
                features[validation_index],
                features[test_index],
                fold,
                c,
            )
            oof_predictions.setdefault(method, np.empty(len(y)))[validation] = val_score
            test_predictions.setdefault(method, []).append(test_score)

    scaled_embeddings = {}
    for embedding_name, features in embeddings.items():
        scaler = StandardScaler().fit(features[train_index])
        scaled_embeddings[embedding_name] = (
            scaler.transform(features[train_index]),
            scaler.transform(features[validation_index]),
            scaler.transform(features[test_index]),
        )
        for c in [0.02, 0.05, 0.2]:
            method = f"logistic_{embedding_name}_c{c:g}"
            val_score, test_score = fit_logistic(
                method,
                scaled_embeddings[embedding_name][0],
                y[train],
                scaled_embeddings[embedding_name][1],
                scaled_embeddings[embedding_name][2],
                fold,
                c,
            )
            oof_predictions.setdefault(method, np.empty(len(y)))[validation] = val_score
            test_predictions.setdefault(method, []).append(test_score)

    for fp_name, fp_features in fingerprints.items():
        for embedding_name, (emb_train, emb_val, emb_test) in scaled_embeddings.items():
            method = f"logistic_{fp_name}_{embedding_name}_c0.05"
            val_score, test_score = fit_logistic(
                method,
                np.column_stack([fp_features[train_index], 0.5 * emb_train]),
                y[train],
                np.column_stack([fp_features[validation_index], 0.5 * emb_val]),
                np.column_stack([fp_features[test_index], 0.5 * emb_test]),
                fold,
            )
            oof_predictions.setdefault(method, np.empty(len(y)))[validation] = val_score
            test_predictions.setdefault(method, []).append(test_score)
    print(f"fold {fold} complete", flush=True)

current_oof = (
    read("current_oof")
    .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == "all_plain_pfn_graph_mean"))
    .sort("compound_id")
)
current_test = (
    read("current_test")
    .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == "all_plain_pfn_graph_mean"))
    .group_by("compound_id")
    .agg(pl.col("prediction").mean())
    .join(test, on="compound_id", validate="1:1", maintain_order="right")
)
assert current_oof["compound_id"].to_list() == data["compound_id"].to_list()
assert current_test["compound_id"].to_list() == test["compound_id"].to_list()

ranked_oof = {}
raw = {"current": current_oof["prediction"].to_numpy(), **oof_predictions}
for method, prediction in raw.items():
    ranked = np.empty(len(y))
    for fold in range(5):
        mask = folds == fold
        ranked[mask] = unit_rank(prediction[mask])
    ranked_oof[method] = ranked
ranked_test = {"current": unit_rank(current_test["prediction"].to_numpy())}
for method, predictions in test_predictions.items():
    fold_ranks = np.stack([unit_rank(prediction) for prediction in predictions])
    ranked_test[method] = unit_rank(fold_ranks.mean(axis=0))
for method in oof_predictions:
    blend_name = f"current2_{method}"
    ranked_oof[blend_name] = (2 * ranked_oof["current"] + ranked_oof[method]) / 3
    ranked_test[blend_name] = (2 * ranked_test["current"] + ranked_test[method]) / 3

rates = cfg["rates"]
screen_rows = []
for method, prediction in ranked_oof.items():
    for rate_name, rate in rates.items():
        calls = fixed_rate_calls(prediction, folds, rate)
        screen_rows.append(
            {
                "method": method,
                "rate_policy": rate_name,
                "rate": rate,
                "source_MCC": float(matthews_corrcoef(y, calls)),
                "expected_target_MCC": expected_mcc(y, calls, target_prior),
                "AUROC": float(roc_auc_score(y, prediction)),
                "AP": float(average_precision_score(y, prediction)),
            }
        )
screen = pl.DataFrame(screen_rows).sort(
    "rate_policy", "expected_target_MCC", descending=[False, True]
)
screen.write_csv(out / "screen.csv")

selection_rows = []
deployment = {}
methods = list(ranked_oof)
for rate_name, rate in rates.items():
    for top_k in [1, 3, 5, 10]:
        crossfit_prediction = np.empty(len(y))
        fold_selections = []
        for heldout in range(5):
            selection_mask = folds != heldout
            selection_folds = folds[selection_mask]
            scored = []
            for method in methods:
                calls = fixed_rate_calls(ranked_oof[method][selection_mask], selection_folds, rate)
                scored.append(
                    (
                        expected_mcc(y[selection_mask], calls, target_prior),
                        method,
                    )
                )
            chosen = [method for _, method in sorted(scored, reverse=True)[:top_k]]
            fold_selections.append(chosen)
            mask = folds == heldout
            crossfit_prediction[mask] = np.mean(
                [ranked_oof[method][mask] for method in chosen], axis=0
            )
        crossfit_calls = fixed_rate_calls(crossfit_prediction, folds, rate)
        full_scores = []
        for method in methods:
            calls = fixed_rate_calls(ranked_oof[method], folds, rate)
            full_scores.append((expected_mcc(y, calls, target_prior), method))
        deployment_methods = [method for _, method in sorted(full_scores, reverse=True)[:top_k]]
        deployment_name = f"{rate_name}_top{top_k}"
        deployment[deployment_name] = {
            "methods": deployment_methods,
            "prediction": np.mean([ranked_test[method] for method in deployment_methods], axis=0),
            "rate": rate,
        }
        selection_rows.append(
            {
                "variant": deployment_name,
                "rate": rate,
                "top_k": top_k,
                "crossfit_source_MCC": float(matthews_corrcoef(y, crossfit_calls)),
                "crossfit_expected_target_MCC": expected_mcc(y, crossfit_calls, target_prior),
                "crossfit_AUROC": float(roc_auc_score(y, crossfit_prediction)),
                "deployment_methods": json.dumps(deployment_methods),
                "fold_methods": json.dumps(fold_selections),
            }
        )
selection = pl.DataFrame(selection_rows).sort("crossfit_expected_target_MCC", descending=True)
selection.write_csv(out / "crossfit-selection.csv")

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
test_rows = []
for method, prediction in ranked_test.items():
    test_rows.append(
        test.with_columns(
            pl.lit(endpoint).alias("endpoint"),
            pl.lit(method).alias("method"),
            pl.Series("prediction", prediction),
        )
    )
for variant, values in deployment.items():
    calls = np.zeros(len(test), dtype=bool)
    count = int(round(values["rate"] * len(test)))
    calls[np.argsort(-values["prediction"], kind="stable")[:count]] = True
    test_rows.append(
        test.with_columns(
            pl.lit(endpoint).alias("endpoint"),
            pl.lit(variant).alias("method"),
            pl.Series("prediction", values["prediction"]),
            pl.Series("binary", calls),
        )
    )
pl.concat(test_rows, how="diagonal").write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "endpoint": endpoint,
            "fingerprints": list(fingerprints),
            "embeddings": list(embeddings),
            "base_models": len(oof_predictions),
            "rank_candidates": len(ranked_oof),
            "selection": "For each heldout outer fold, rank methods on the other four OOF folds.",
            "target_prior": target_prior,
        },
        indent=2,
    )
    + "\n"
)
print(selection.head(12).write_csv(), flush=True)
