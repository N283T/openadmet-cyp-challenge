"""Build leakage-resistant classical PXR predictions for CYP compounds."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GroupKFold

from cyp_challenge.preprocessing import standardize

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
RDLogger.DisableLog("rdApp.*")


def normalize(smiles: str) -> tuple[str | None, str | None]:
    result = standardize(smiles)
    normalized = result["standardized_smiles"]
    molecule = Chem.MolFromSmiles(normalized) if normalized else None
    key = Chem.MolToInchiKey(molecule).split("-")[0] if molecule else None
    return normalized, key


def prepare_source(table: pl.DataFrame, smiles_col: str, value_col: str) -> pl.DataFrame:
    rows = []
    for smiles, value in table.select(smiles_col, value_col).iter_rows():
        normalized, key = normalize(smiles)
        if key is not None and value is not None and np.isfinite(value):
            rows.append({"connectivity": key, "smiles": normalized, "pec50": float(value)})
    return pl.DataFrame(rows)


def scaffold(smiles: str, key: str) -> str:
    value = MurckoScaffold.MurckoScaffoldSmiles(smiles=smiles)
    return value or f"acyclic:{key}"


def fingerprints(smiles: list[str], radius: int, bits: int) -> np.ndarray:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=bits)
    matrix = np.zeros((len(smiles), bits), dtype=np.float32)
    for index, value in enumerate(smiles):
        fp = generator.GetFingerprint(Chem.MolFromSmiles(value))
        DataStructs.ConvertToNumpyArray(fp, matrix[index])
    return matrix


def models(seed: int) -> dict[str, object]:
    return {
        "ridge": Ridge(alpha=10.0),
        "extra_trees": ExtraTreesRegressor(
            n_estimators=500,
            min_samples_leaf=2,
            max_features=0.5,
            n_jobs=-1,
            random_state=seed,
        ),
        "random_forest": RandomForestRegressor(
            n_estimators=400,
            min_samples_leaf=2,
            max_features=0.5,
            n_jobs=-1,
            random_state=seed,
        ),
        "lgbm": lgb.LGBMRegressor(
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
    }


radius = int(cfg["fingerprint_radius"])
bits = int(cfg["fingerprint_bits"])
seed = int(cfg["random_seed"])
compounds = pl.read_parquet(inp / cfg["compounds"]).select("compound_id", "standardized_smiles")
challenge_keys = []
for smiles in compounds["standardized_smiles"]:
    molecule = Chem.MolFromSmiles(smiles)
    challenge_keys.append(Chem.MolToInchiKey(molecule).split("-")[0])
challenge_key_set = set(challenge_keys)

external = pl.concat(
    [
        prepare_source(pl.read_parquet(inp / cfg["pxr_train_activity"]), "smiles", "pec50"),
        prepare_source(pl.read_parquet(inp / cfg["pxr_test_activity"]), "SMILES", "pEC50"),
    ]
)
external = (
    external.filter(~pl.col("connectivity").is_in(challenge_key_set))
    .group_by("connectivity")
    .agg(pl.col("smiles").first(), pl.col("pec50").mean())
    .sort("connectivity")
)
external = external.with_columns(
    pl.struct("smiles", "connectivity")
    .map_elements(lambda row: scaffold(row["smiles"], row["connectivity"]), return_dtype=pl.String)
    .alias("scaffold")
)

x_train = fingerprints(external["smiles"].to_list(), radius, bits)
x_cyp = fingerprints(compounds["standardized_smiles"].to_list(), radius, bits)
y = external["pec50"].to_numpy()
groups = external["scaffold"].to_numpy()
splitter = GroupKFold(n_splits=int(cfg["cv_folds"]))
oof_columns: dict[str, np.ndarray] = {}
metrics = []
for name in models(seed):
    oof = np.full(len(y), np.nan)
    for fold, (train_index, valid_index) in enumerate(splitter.split(x_train, y, groups)):
        model = models(seed + fold)[name]
        model.fit(x_train[train_index], y[train_index])
        oof[valid_index] = model.predict(x_train[valid_index])
    oof_columns[name] = oof
    metrics.append(
        {
            "model": name,
            "rmse": float(mean_squared_error(y, oof) ** 0.5),
            "mae": float(mean_absolute_error(y, oof)),
            "spearman": float(spearmanr(y, oof).statistic),
        }
    )

predictions = {}
for name, model in models(seed).items():
    model.fit(x_train, y)
    predictions[name] = model.predict(x_cyp)
predictions["mean_tree"] = np.mean(
    np.column_stack(
        [predictions["extra_trees"], predictions["random_forest"], predictions["lgbm"]]
    ),
    axis=1,
)

feature_table = compounds.select("compound_id").with_columns(
    [pl.Series(f"pxr_pec50_{name}", values) for name, values in predictions.items()]
)
feature_table.write_parquet(out / "pxr-classical-features.parquet")
oof_table = external.select("connectivity", "pec50", "scaffold").with_columns(
    [pl.Series(f"oof_{name}", values) for name, values in oof_columns.items()]
)
oof_table.write_parquet(out / "pxr-external-oof.parquet")
pl.DataFrame(metrics).sort("rmse").write_csv(out / "cv-metrics.csv")
summary = {
    "external_rows_after_cyp_exclusion": external.height,
    "excluded_cyp_connectivities": len(challenge_key_set),
    "external_scaffolds": external["scaffold"].n_unique(),
    "feature_columns": [column for column in feature_table.columns if column != "compound_id"],
    "cv_metrics": sorted(metrics, key=lambda row: row["rmse"]),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
