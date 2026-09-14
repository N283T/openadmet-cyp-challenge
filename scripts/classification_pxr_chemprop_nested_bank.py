"""Build outer/inner/test probabilities for strict PXR Chemprop readouts."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
RDLogger.DisableLog("rdApp.*")


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def fit_predict(
    embedding: np.ndarray,
    morgan: np.ndarray,
    train_index: np.ndarray,
    y: np.ndarray,
    queries: list[np.ndarray],
    seed: int,
) -> dict[str, list[np.ndarray]]:
    x_train = embedding[train_index]
    query_embedding = [embedding[index] for index in queries]
    models = {
        "strict_cp_lgbm": lgb.LGBMClassifier(
            n_estimators=500,
            learning_rate=0.03,
            num_leaves=31,
            min_child_samples=20,
            subsample=0.8,
            colsample_bytree=0.6,
            reg_lambda=1.0,
            verbosity=-1,
            n_jobs=-1,
            random_state=seed,
        ),
        "strict_cp_extra_trees": ExtraTreesClassifier(
            n_estimators=400,
            min_samples_leaf=2,
            max_features="sqrt",
            n_jobs=-1,
            random_state=seed,
        ),
    }
    result = {}
    for name, model in models.items():
        model.fit(x_train, y)
        result[name] = [model.predict_proba(query)[:, 1] for query in query_embedding]

    scaler = StandardScaler().fit(x_train)
    svm = SVC(C=1.0, probability=True, random_state=seed).fit(scaler.transform(x_train), y)
    result["strict_cp_rbf_svm"] = [
        svm.predict_proba(scaler.transform(query))[:, 1] for query in query_embedding
    ]

    combined_train = np.column_stack([morgan[train_index], x_train])
    combined_queries = [
        np.column_stack([morgan[index], query])
        for index, query in zip(queries, query_embedding, strict=True)
    ]
    combined = lgb.LGBMClassifier(
        n_estimators=500,
        learning_rate=0.03,
        num_leaves=31,
        min_child_samples=20,
        subsample=0.8,
        colsample_bytree=0.6,
        reg_lambda=1.0,
        verbosity=-1,
        n_jobs=-1,
        random_state=seed,
    ).fit(combined_train, y)
    result["strict_morgan_cp_lgbm"] = [
        combined.predict_proba(query)[:, 1] for query in combined_queries
    ]
    result["strict_cp_threeway"] = [
        np.mean(
            np.column_stack(
                [
                    result["strict_morgan_cp_lgbm"][query_index],
                    result["strict_cp_extra_trees"][query_index],
                    result["strict_cp_rbf_svm"][query_index],
                ]
            ),
            axis=1,
        )
        for query_index in range(len(queries))
    ]
    return result


compounds = read("compounds").sort("compound_id")
all_ids = compounds.select("compound_id")
id_to_index = {value: row for row, value in enumerate(all_ids["compound_id"])}
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
test_index = np.asarray([id_to_index[value] for value in test["compound_id"]])
feature_table = read("features")
columns = [column for column in feature_table.columns if column.startswith("pxr_cp_emb_")]
embedding = (
    all_ids.join(feature_table.select("compound_id", *columns), on="compound_id", validate="1:1")
    .select(columns)
    .to_numpy()
    .astype(np.float32)
)
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in compounds["standardized_smiles"]
    ]
).astype(np.float32)
assignments = read("assignments").filter(pl.col("strategy") == cfg["strategy"])
labels = read("labels")

outer_rows = []
inner_rows = []
test_rows = []
for endpoint in cfg["endpoints"]:
    data = (
        labels.filter(pl.col("endpoint") == endpoint)
        .select("compound_id", "value")
        .join(assignments.select("compound_id", "fold"), on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    data_index = np.asarray([id_to_index[value] for value in data["compound_id"]])
    y = data["value"].to_numpy().astype(int)
    folds = data["fold"].to_numpy().astype(int)
    for outer_fold in np.unique(folds):
        train_rows = np.flatnonzero(folds != outer_fold)
        heldout_rows = np.flatnonzero(folds == outer_fold)
        prediction = fit_predict(
            embedding,
            morgan,
            data_index[train_rows],
            y[train_rows],
            [data_index[heldout_rows], test_index],
            cfg["random_seed"] + int(outer_fold),
        )
        for method, (heldout_prediction, fold_test_prediction) in prediction.items():
            for row, probability in zip(heldout_rows, heldout_prediction, strict=True):
                outer_rows.append(
                    {
                        "compound_id": data["compound_id"][int(row)],
                        "endpoint": endpoint,
                        "fold": int(outer_fold),
                        "y_true": int(y[row]),
                        "probability": float(probability),
                        "method": method,
                        "split_id": cfg["split_id"],
                    }
                )
            for row, probability in enumerate(fold_test_prediction):
                test_rows.append(
                    {
                        "compound_id": test["compound_id"][row],
                        "endpoint": endpoint,
                        "method": method,
                        "fold": int(outer_fold),
                        "split_id": cfg["split_id"],
                        "probability": float(probability),
                    }
                )

    for inner_fold in np.unique(folds):
        outer_fold = (int(inner_fold) - 1) % 5
        train_rows = np.flatnonzero((folds != inner_fold) & (folds != outer_fold))
        heldout_rows = np.flatnonzero(folds == inner_fold)
        prediction = fit_predict(
            embedding,
            morgan,
            data_index[train_rows],
            y[train_rows],
            [data_index[heldout_rows]],
            cfg["random_seed"] + 100 + int(inner_fold),
        )
        for method, (heldout_prediction,) in prediction.items():
            for row, probability in zip(heldout_rows, heldout_prediction, strict=True):
                inner_rows.append(
                    {
                        "compound_id": data["compound_id"][int(row)],
                        "endpoint": endpoint,
                        "fold": int(inner_fold),
                        "y_true": int(y[row]),
                        "probability": float(probability),
                        "method": method,
                        "split_id": cfg["split_id"],
                        "outer_fold": outer_fold,
                        "inner_fold": int(inner_fold),
                    }
                )
        print(f"{endpoint} inner fold {inner_fold} complete", flush=True)

pl.DataFrame(outer_rows).write_parquet(out / "outer-oof.parquet")
pl.DataFrame(inner_rows).write_parquet(out / "inner-probabilities.parquet")
pl.DataFrame(test_rows).write_parquet(out / "fold-test-probabilities.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "methods": sorted(set(row["method"] for row in outer_rows)),
            "outer_rows": len(outer_rows),
            "inner_rows": len(inner_rows),
            "test_rows": len(test_rows),
            "inner_rotation": "outer_fold = (inner_fold - 1) mod 5",
        },
        indent=2,
    )
    + "\n"
)
