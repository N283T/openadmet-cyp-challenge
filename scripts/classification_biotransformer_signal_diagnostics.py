"""Diagnose target association and importance of BioTransformer-derived features."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.metrics import average_precision_score, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def activity_block(name: str) -> str:
    if name.endswith("_parent"):
        return "activity_parent"
    scope = next(scope for scope in ("all", "d6", "a4") if f"_{scope}_" in name)
    kind = "delta" if name.endswith("_delta") else "absolute"
    return f"activity_{scope}_{kind}"


def association_rows(
    matrix: np.ndarray,
    names: list[str],
    y: np.ndarray,
    endpoint: str,
    strategy: str,
    source: str,
) -> list[dict[str, object]]:
    rows = []
    for column, name in enumerate(names):
        values = matrix[:, column]
        if np.std(values) == 0:
            pearson = spearman = 0.0
        else:
            pearson = float(pearsonr(values, y).statistic)
            spearman = float(spearmanr(values, y).statistic)
        rows.append(
            {
                "endpoint": endpoint,
                "strategy": strategy,
                "source": source,
                "feature": name,
                "pearson": pearson,
                "spearman": spearman,
                "abs_spearman": abs(spearman),
                "positive_mean": float(values[y == 1].mean()),
                "negative_mean": float(values[y == 0].mean()),
            }
        )
    return rows


def fit_model(x: np.ndarray, y: np.ndarray, fold: int) -> lgb.LGBMClassifier:
    model = lgb.LGBMClassifier(
        **cfg["lgbm"], verbosity=-1, n_jobs=1, random_state=cfg["seed"] + fold
    )
    model.fit(x, y)
    return model


def permutation_rows(
    model: lgb.LGBMClassifier,
    x: np.ndarray,
    y: np.ndarray,
    blocks: dict[str, np.ndarray],
    endpoint: str,
    strategy: str,
    method: str,
    fold: int,
) -> list[dict[str, object]]:
    base = model.predict_proba(x)[:, 1]
    base_auroc = roc_auc_score(y, base)
    base_ap = average_precision_score(y, base)
    rng = np.random.default_rng(cfg["seed"] + 1009 * fold)
    rows = []
    for block, columns in blocks.items():
        auroc_drops = []
        ap_drops = []
        for _ in range(cfg["permutation_repeats"]):
            shuffled = x.copy()
            order = rng.permutation(len(x))
            shuffled[:, columns] = x[order][:, columns]
            score = model.predict_proba(shuffled)[:, 1]
            auroc_drops.append(base_auroc - roc_auc_score(y, score))
            ap_drops.append(base_ap - average_precision_score(y, score))
        rows.append(
            {
                "endpoint": endpoint,
                "strategy": strategy,
                "method": method,
                "fold": fold,
                "block": block,
                "AUROC_drop": float(np.mean(auroc_drops)),
                "AP_drop": float(np.mean(ap_drops)),
            }
        )
    return rows


compounds = read("compounds").sort("compound_id")
labels = read("labels")
assignments = read("assignments")
activity = read("activity_features")
embeddings = read("embeddings").sort("compound_id")
embedding_scalars = read("embedding_scalars").sort("compound_id")
if not embeddings.select("compound_id").equals(compounds.select("compound_id")):
    raise ValueError("Embedding identities differ from compounds")
if not embedding_scalars.select("compound_id").equals(compounds.select("compound_id")):
    raise ValueError("Embedding scalar identities differ from compounds")

id_to_index = {compound_id: index for index, compound_id in enumerate(compounds["compound_id"])}
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in compounds["standardized_smiles"]
    ]
).astype(np.float32)


def embedding_block(scope: str) -> np.ndarray:
    columns = sorted(
        column for column in embeddings.columns if column.startswith(f"chemeleon_{scope}_")
    )
    return embeddings.select(columns).to_numpy().astype(np.float32)


parent = embedding_block("parent")
scope_embeddings = {"d6": embedding_block("d6"), "a4": embedding_block("a4")}
scalar_names = [column for column in embedding_scalars.columns if column != "compound_id"]
scalar_matrix = embedding_scalars.select(scalar_names).to_numpy().astype(np.float32)

association = []
feature_importance = []
permutation = []
for endpoint in cfg["endpoints"]:
    enzyme_scope = {"CYP2D6": "d6", "CYP3A4": "a4"}[endpoint.split("_", 1)[0]]
    product_delta = scope_embeddings[enzyme_scope] - parent
    endpoint_labels = labels.filter(pl.col("endpoint") == endpoint).select(
        "compound_id", "value"
    )
    for strategy in cfg["strategies"]:
        fold_table = assignments.filter(pl.col("strategy") == strategy).select(
            "compound_id", "fold"
        )
        data = endpoint_labels.join(fold_table, on="compound_id", validate="1:1").sort(
            "compound_id"
        )
        indices = np.asarray([id_to_index[compound_id] for compound_id in data["compound_id"]])
        folds = data["fold"].to_numpy().astype(int)
        y = data["value"].to_numpy().astype(int)
        activity_names: list[str] | None = None
        oof_activity: np.ndarray | None = None
        for fold in range(5):
            train_rows = np.flatnonzero(folds != fold)
            heldout_rows = np.flatnonzero(folds == fold)
            fold_features = (
                activity.filter(
                    (pl.col("strategy") == strategy) & (pl.col("model_fold") == fold)
                )
                .drop("strategy", "model_fold")
                .sort("compound_id")
            )
            if not fold_features.select("compound_id").equals(
                compounds.select("compound_id")
            ):
                raise ValueError("Activity feature grid does not match compounds")
            names = [column for column in fold_features.columns if column != "compound_id"]
            if activity_names is None:
                activity_names = names
                oof_activity = np.empty((data.height, len(names)), dtype=np.float32)
            if names != activity_names or oof_activity is None:
                raise ValueError("Activity feature schema changed across folds")
            activity_matrix = fold_features.select(names).to_numpy().astype(np.float32)
            oof_activity[heldout_rows] = activity_matrix[indices[heldout_rows]]

            activity_x = np.column_stack([morgan[indices], activity_matrix[indices]])
            activity_model = fit_model(activity_x[train_rows], y[train_rows], fold)
            gain = activity_model.booster_.feature_importance(importance_type="gain")
            split = activity_model.booster_.feature_importance(importance_type="split")
            gain_total = max(float(gain.sum()), 1.0)
            split_total = max(float(split.sum()), 1.0)
            feature_names = [f"morgan_{index}" for index in range(2048)] + names
            for feature, gain_value, split_value in zip(
                feature_names, gain, split, strict=True
            ):
                feature_importance.append(
                    {
                        "endpoint": endpoint,
                        "strategy": strategy,
                        "method": "morgan_activity_all_lgbm",
                        "fold": fold,
                        "feature": feature,
                        "block": (
                            "morgan"
                            if feature.startswith("morgan_")
                            else activity_block(feature)
                        ),
                        "gain_fraction": float(gain_value / gain_total),
                        "split_fraction": float(split_value / split_total),
                    }
                )
            block_columns: dict[str, list[int]] = {"morgan": list(range(2048))}
            for column, name in enumerate(names, start=2048):
                block_columns.setdefault(activity_block(name), []).append(column)
            permutation.extend(
                permutation_rows(
                    activity_model,
                    activity_x[heldout_rows],
                    y[heldout_rows],
                    {key: np.asarray(value) for key, value in block_columns.items()},
                    endpoint,
                    strategy,
                    "morgan_activity_all_lgbm",
                    fold,
                )
            )

            for method, raw in (
                ("morgan_parent_pca_lgbm", parent),
                ("morgan_product_delta_endpoint_pca_lgbm", product_delta),
            ):
                pca = PCA(
                    n_components=cfg["pca_components"],
                    svd_solver="randomized",
                    iterated_power=3,
                    random_state=cfg["seed"] + fold,
                )
                pca.fit(raw[indices[train_rows]])
                reduced = pca.transform(raw[indices]).astype(np.float32)
                x = np.column_stack([morgan[indices], reduced])
                model = fit_model(x[train_rows], y[train_rows], fold)
                gain = model.booster_.feature_importance(importance_type="gain")
                split = model.booster_.feature_importance(importance_type="split")
                gain_total = max(float(gain.sum()), 1.0)
                split_total = max(float(split.sum()), 1.0)
                for block, columns in {
                    "morgan": np.arange(2048),
                    "embedding": np.arange(2048, 2048 + cfg["pca_components"]),
                }.items():
                    feature_importance.append(
                        {
                            "endpoint": endpoint,
                            "strategy": strategy,
                            "method": method,
                            "fold": fold,
                            "feature": f"{block}_aggregate",
                            "block": block,
                            "gain_fraction": float(gain[columns].sum() / gain_total),
                            "split_fraction": float(split[columns].sum() / split_total),
                        }
                    )
                permutation.extend(
                    permutation_rows(
                        model,
                        x[heldout_rows],
                        y[heldout_rows],
                        {
                            "morgan": np.arange(2048),
                            "embedding": np.arange(
                                2048, 2048 + cfg["pca_components"]
                            ),
                        },
                        endpoint,
                        strategy,
                        method,
                        fold,
                    )
                )
        if activity_names is None or oof_activity is None:
            raise ValueError("No activity features were routed")
        association.extend(
            association_rows(
                oof_activity, activity_names, y, endpoint, strategy, "activity_oof"
            )
        )
        association.extend(
            association_rows(
                scalar_matrix[indices],
                scalar_names,
                y,
                endpoint,
                strategy,
                "embedding_scalar",
            )
        )

association_table = pl.DataFrame(association).sort(
    "endpoint", "strategy", "abs_spearman", descending=[False, False, True]
)
importance_table = pl.DataFrame(feature_importance)
importance_by_fold = importance_table.group_by(
    "endpoint", "strategy", "method", "fold", "block"
).agg(
    pl.col("gain_fraction").sum().alias("gain_fraction"),
    pl.col("split_fraction").sum().alias("split_fraction"),
)
importance_summary = (
    importance_by_fold.group_by("endpoint", "strategy", "method", "block")
    .agg(
        pl.col("gain_fraction").mean().alias("mean_gain_fraction"),
        pl.col("split_fraction").mean().alias("mean_split_fraction"),
    )
    .sort("endpoint", "strategy", "method", "mean_gain_fraction", descending=True)
)
permutation_table = pl.DataFrame(permutation)
permutation_summary = (
    permutation_table.group_by("endpoint", "strategy", "method", "block")
    .agg(
        pl.col("AUROC_drop").mean().alias("mean_AUROC_drop"),
        pl.col("AP_drop").mean().alias("mean_AP_drop"),
    )
    .sort("endpoint", "strategy", "method", "mean_AUROC_drop", descending=True)
)
association_table.write_csv(out / "target-associations.csv")
importance_table.write_parquet(out / "feature-importance-by-fold.parquet")
importance_summary.write_csv(out / "importance-by-block.csv")
permutation_table.write_csv(out / "permutation-by-fold.csv")
permutation_summary.write_csv(out / "permutation-by-block.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "selection": cfg["selection"],
            "association_rows": association_table.height,
            "importance_rows": importance_table.height,
            "permutation_rows": permutation_table.height,
            "notes": [
                "Activity associations use transforms excluding each row's outer fold.",
                "Gain/split fractions are normalized per model before fold averaging.",
                "Permutation importance is computed only on each held-out fold.",
                "Embedding PCA axes are not compared across folds; importance is block-aggregated.",
            ],
        },
        indent=2,
    )
    + "\n"
)
print(association_table.group_by("endpoint", "strategy").head(8), flush=True)
print(importance_summary, flush=True)
print(permutation_summary, flush=True)
