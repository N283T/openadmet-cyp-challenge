"""Evaluate PXR and other compact predicted phenotypes for CYP3A4 TDI."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from sklearn.preprocessing import StandardScaler

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
RDLogger.DisableLog("rdApp.*")


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def align(all_ids: pl.DataFrame, table: pl.DataFrame, columns: list[str]) -> np.ndarray:
    joined = all_ids.join(
        table.select("compound_id", *columns),
        on="compound_id",
        how="left",
        validate="1:1",
        maintain_order="left",
    )
    values = joined.select(columns).to_numpy().astype(np.float32)
    if not np.isfinite(values).all():
        raise ValueError(f"Incomplete auxiliary features: {columns}")
    return values


def fixed_rate_calls(scores: np.ndarray, folds: np.ndarray, rate: float) -> np.ndarray:
    calls = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = max(1, int(round(rate * len(positions))))
        calls[positions[np.argsort(-scores[positions], kind="stable")[:count]]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray, prevalence: float) -> float:
    positive = y == 1
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    tp = prevalence * tpr
    fn = prevalence * (1 - tpr)
    fp = (1 - prevalence) * fpr
    tn = (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def estimator(kind: str, seed: int, positive_weight: float):
    if kind == "lr":
        return LogisticRegression(C=0.05, solver="liblinear", max_iter=4000, random_state=seed)
    if kind == "lr_balanced":
        return LogisticRegression(
            C=0.05,
            solver="liblinear",
            max_iter=4000,
            class_weight="balanced",
            random_state=seed,
        )
    parameters = {
        "n_estimators": 500,
        "learning_rate": 0.03,
        "num_leaves": 31,
        "min_child_samples": 20,
        "subsample": 0.8,
        "colsample_bytree": 0.6,
        "reg_lambda": 1.0,
        "verbosity": -1,
        "n_jobs": -1,
        "random_state": seed,
    }
    if kind == "lgbm_balanced":
        parameters["scale_pos_weight"] = positive_weight
    return lgb.LGBMClassifier(**parameters)


compounds = read("compounds").sort("compound_id")
all_ids = compounds.select("compound_id")
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
id_to_index = {value: row for row, value in enumerate(all_ids["compound_id"])}
test_index = np.asarray([id_to_index[value] for value in test["compound_id"]])

generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in compounds["standardized_smiles"]
    ]
).astype(np.float32)
pxr_table = read("pxr_features")
pxr_columns = [column for column in pxr_table.columns if column.startswith("pxr_")]
pxr = align(all_ids, pxr_table, pxr_columns)
admetica_table = read("admetica_features")
admetica_columns = [column for column in admetica_table.columns if column.endswith("_probability")]
admetica = align(all_ids, admetica_table, admetica_columns)
substrate_table = read("substrate_features")
substrate_columns = [
    column for column in substrate_table.columns if column.startswith("substrate_")
]
substrate = align(all_ids, substrate_table, substrate_columns)
pred_bank = np.column_stack([pxr, admetica, substrate])
aux_no_pxr = np.column_stack([admetica, substrate])

feature_views = {
    "pxr": pxr,
    "pred_bank": pred_bank,
    "morgan": morgan,
    "morgan_pxr": np.column_stack([morgan, pxr]),
    "morgan_admetica": np.column_stack([morgan, admetica]),
    "morgan_substrate": np.column_stack([morgan, substrate]),
    "morgan_aux_no_pxr": np.column_stack([morgan, aux_no_pxr]),
    "morgan_pxr_admetica": np.column_stack([morgan, pxr, admetica]),
    "morgan_pxr_substrate": np.column_stack([morgan, pxr, substrate]),
    "morgan_pred_bank": np.column_stack([morgan, pred_bank]),
}
methods = [
    ("pxr_lr", "pxr", "lr"),
    ("pxr_lr_balanced", "pxr", "lr_balanced"),
    ("pred_bank_lr", "pred_bank", "lr"),
    ("morgan_lr", "morgan", "lr"),
    ("morgan_pxr_lr", "morgan_pxr", "lr"),
    ("morgan_pred_bank_lr", "morgan_pred_bank", "lr"),
    ("morgan_lgbm", "morgan", "lgbm"),
    ("morgan_pxr_lgbm", "morgan_pxr", "lgbm"),
    ("morgan_pxr_lgbm_balanced", "morgan_pxr", "lgbm_balanced"),
    ("morgan_admetica_lgbm", "morgan_admetica", "lgbm"),
    ("morgan_substrate_lgbm", "morgan_substrate", "lgbm"),
    ("morgan_aux_no_pxr_lgbm", "morgan_aux_no_pxr", "lgbm"),
    ("morgan_pxr_admetica_lgbm", "morgan_pxr_admetica", "lgbm"),
    ("morgan_pxr_substrate_lgbm", "morgan_pxr_substrate", "lgbm"),
    ("morgan_pred_bank_lgbm", "morgan_pred_bank", "lgbm"),
    ("morgan_pred_bank_lgbm_balanced", "morgan_pred_bank", "lgbm_balanced"),
]
direct_blocks = [(pxr, pxr_columns), (admetica, admetica_columns), (substrate, substrate_columns)]
if "pxr_chemprop_features" in cfg:
    chemprop_table = read("pxr_chemprop_features")
    chemprop_embedding_columns = [
        column for column in chemprop_table.columns if column.startswith("pxr_cp_emb_")
    ]
    chemprop_scalar_columns = [
        column
        for column in chemprop_table.columns
        if column.startswith("pxr_cp_") and column not in chemprop_embedding_columns
    ]
    chemprop_embedding = align(all_ids, chemprop_table, chemprop_embedding_columns)
    chemprop_scalar = align(all_ids, chemprop_table, chemprop_scalar_columns)
    feature_views.update(
        {
            "chemprop_scalar": chemprop_scalar,
            "chemprop_embedding": chemprop_embedding,
            "morgan_chemprop_scalar": np.column_stack([morgan, chemprop_scalar]),
            "morgan_chemprop_embedding": np.column_stack([morgan, chemprop_embedding]),
            "morgan_pxr_admetica_cp_scalar": np.column_stack(
                [morgan, pxr, admetica, chemprop_scalar]
            ),
            "morgan_pxr_admetica_cp_embedding": np.column_stack(
                [morgan, pxr, admetica, chemprop_embedding]
            ),
        }
    )
    methods.extend(
        [
            ("chemprop_scalar_lr", "chemprop_scalar", "lr"),
            ("chemprop_embedding_lr", "chemprop_embedding", "lr"),
            ("chemprop_embedding_lgbm", "chemprop_embedding", "lgbm"),
            ("morgan_chemprop_scalar_lgbm", "morgan_chemprop_scalar", "lgbm"),
            ("morgan_chemprop_embedding_lgbm", "morgan_chemprop_embedding", "lgbm"),
            (
                "morgan_pxr_admetica_cp_scalar_lgbm",
                "morgan_pxr_admetica_cp_scalar",
                "lgbm",
            ),
            (
                "morgan_pxr_admetica_cp_embedding_lgbm",
                "morgan_pxr_admetica_cp_embedding",
                "lgbm",
            ),
        ]
    )
    direct_blocks.append((chemprop_scalar, chemprop_scalar_columns))

labels = read("labels").filter(pl.col("endpoint") == cfg["endpoint"]).select("compound_id", "value")
metric_rows = []
oof_rows = []
test_rows = []
for strategy in cfg["strategies"]:
    fold_table = (
        read("assignments")
        .filter(pl.col("strategy") == strategy)
        .select("compound_id", "fold", "group")
    )
    data = labels.join(fold_table, on="compound_id", validate="1:1").sort("compound_id")
    data_index = np.asarray([id_to_index[value] for value in data["compound_id"]])
    y = data["value"].to_numpy().astype(int)
    folds = data["fold"].to_numpy().astype(int)
    positive_weight = float((y == 0).sum() / (y == 1).sum())
    predictions = {name: np.empty(len(data)) for name, _, _ in methods}
    test_predictions = {name: [] for name, _, _ in methods}
    for fold in sorted(set(folds)):
        train_rows = np.flatnonzero(folds != fold)
        valid_rows = np.flatnonzero(folds == fold)
        for name, view, kind in methods:
            scaler = StandardScaler() if kind.startswith("lr") else None
            train_features = feature_views[view][data_index[train_rows]]
            valid_features = feature_views[view][data_index[valid_rows]]
            fold_test_features = feature_views[view][test_index]
            if scaler is not None:
                train_features = scaler.fit_transform(train_features)
                valid_features = scaler.transform(valid_features)
                fold_test_features = scaler.transform(fold_test_features)
            model = estimator(kind, cfg["random_seed"] + fold, positive_weight)
            model.fit(train_features, y[train_rows])
            predictions[name][valid_rows] = model.predict_proba(valid_features)[:, 1]
            test_predictions[name].append(model.predict_proba(fold_test_features)[:, 1])
        print(f"{strategy} fold {fold} complete", flush=True)

    direct_all = np.column_stack([block for block, _ in direct_blocks])
    direct_columns = [column for _, columns in direct_blocks for column in columns]
    direct_scores = {
        f"direct_{column}": direct_all[data_index, index]
        for index, column in enumerate(direct_columns)
    }
    predictions.update(direct_scores)
    for name, values in predictions.items():
        calls = fixed_rate_calls(values, folds, cfg["call_rate"])
        metric_rows.append(
            {
                "strategy": strategy,
                "method": name,
                "AUROC": float(roc_auc_score(y, values)),
                "AP": float(average_precision_score(y, values)),
                "MCC_observed": float(matthews_corrcoef(y, calls)),
                "MCC_target_prior": expected_mcc(y, calls, cfg["target_prior"]),
                "calls": int(calls.sum()),
            }
        )
        for row, prediction in enumerate(values):
            oof_rows.append(
                {
                    "compound_id": data["compound_id"][row],
                    "endpoint": cfg["endpoint"],
                    "strategy": strategy,
                    "fold": int(folds[row]),
                    "method": name,
                    "prediction": float(prediction),
                    "y_true": int(y[row]),
                }
            )
        if name in test_predictions:
            averaged_test = np.mean(test_predictions[name], axis=0)
        else:
            averaged_test = direct_all[
                test_index, direct_columns.index(name.removeprefix("direct_"))
            ]
        for row, prediction in enumerate(averaged_test):
            test_rows.append(
                {
                    "compound_id": test["compound_id"][row],
                    "endpoint": cfg["endpoint"],
                    "strategy": strategy,
                    "method": name,
                    "prediction": float(prediction),
                    "fold": -1,
                }
            )

metrics = pl.DataFrame(metric_rows).sort("strategy", "MCC_target_prior", descending=[False, True])
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(oof_rows).write_parquet(out / "oof.parquet")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "endpoint": cfg["endpoint"],
            "pxr_columns": pxr_columns,
            "admetica_columns": admetica_columns,
            "substrate_columns": substrate_columns,
            "methods": [name for name, _, _ in methods],
            "note": "Prespecified diagnostic screen; no submission candidate selected.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
