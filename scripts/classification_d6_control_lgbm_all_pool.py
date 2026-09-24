"""Extend frozen September 15 control LightGBM fits to full UMAP folds."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
comp = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
ids = comp.select("compound_id")
assigned = (
    pl.read_parquet(inp / cfg["assignments"])
    .filter(pl.col("strategy") == cfg["split_id"])
    .select("compound_id", "fold", "group")
)
labels = (
    pl.read_parquet(inp / cfg["labels"])
    .filter(pl.col("endpoint") == cfg["endpoint"])
    .select("compound_id", pl.col("value").cast(pl.Int64).alias("y_true"))
)
pool = assigned.join(labels, on="compound_id", how="left", validate="1:1").sort("compound_id")
if pool.height != 6146 or pool["compound_id"].n_unique() != 6146:
    raise ValueError("Incomplete control molecular pool")
lookup = {value: i for i, value in enumerate(comp["compound_id"])}
pool_indices = np.asarray([lookup[value] for value in pool["compound_id"]])
train = labels.join(assigned, on="compound_id", validate="1:1").sort("compound_id")
train_indices = np.asarray([lookup[value] for value in train["compound_id"]])
train_folds = train["fold"].to_numpy()
train_y = train["y_true"].to_numpy().astype(int)


def matrix(table: pl.DataFrame, prefix: str | None = None) -> np.ndarray:
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


mordred = pl.concat(
    [
        pl.read_parquet(inp / f"{cfg['mordred_prefix']}{i:05d}.parquet")
        for i in range(cfg["mordred_files"])
    ]
)
xm = matrix(mordred)
xc = matrix(pl.read_parquet(inp / cfg["class_emb"]), "emb_")
xmix = matrix(pl.read_parquet(inp / cfg["mixed_emb"]), "emb_")
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
xf = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in comp["standardized_smiles"]
    ]
).astype(np.float32)
views = {
    "morgan_plain": xf,
    "mordred_plain": xm,
    "fp_mordred_plain": np.column_stack([xf, xm]),
    "external_class_plain": xc,
    "external_mixed_plain": xmix,
}
parts = []
for method in cfg["methods"]:
    values = views[method]
    for fold in range(5):
        held = pool.filter(pl.col("fold") == fold)
        indices = pool_indices[pool["fold"].to_numpy() == fold]
        alias = f"{cfg['model_prefix']}{method}-fold{fold}.txt"
        if cfg.get("refit", False):
            model = lgb.LGBMClassifier(
                objective="binary",
                n_estimators=400,
                num_leaves=15,
                min_child_samples=20,
                learning_rate=0.03,
                reg_lambda=2,
                n_jobs=3,
                random_state=42 + fold,
                verbosity=-1,
                deterministic=True,
                force_col_wise=True,
            )
            valid = train_folds != fold
            model.fit(values[train_indices[valid]], train_y[valid])
            score = model.predict_proba(values[indices])[:, 1]
        else:
            model = lgb.Booster(model_file=str(inp / alias))
            score = model.predict(values[indices])
        if not np.isfinite(score).all():
            raise ValueError(f"Nonfinite control score: {method} {fold}")
        parts.append(
            held.with_columns(
                pl.lit(cfg["endpoint"]).alias("endpoint"),
                pl.lit(cfg["split_id"]).alias("split_id"),
                pl.lit(method).alias("model"),
                pl.Series("score", score),
            )
        )
    print(f"{method} complete", flush=True)
pl.concat(parts).sort("model", "compound_id").write_parquet(out / "all-pool-oof.parquet")
