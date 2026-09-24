"""Extend frozen Morgan tail models to every UMAP held-fold molecule."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import ComplementNB
from sklearn.svm import SVC

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
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
    raise ValueError("Incomplete fixed-model pool")
lookup = {value: i for i, value in enumerate(compounds["compound_id"])}
train = labels.join(assigned, on="compound_id", validate="1:1").sort("compound_id")
train_indices = np.asarray([lookup[value] for value in train["compound_id"]])
train_folds = train["fold"].to_numpy()
train_y = train["y_true"].to_numpy().astype(int)
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
features = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in compounds["standardized_smiles"]
    ]
).astype(np.float32)
parts = []
for method in cfg["methods"]:
    for fold in range(5):
        held = pool.filter(pl.col("fold") == fold)
        indices = np.asarray([lookup[value] for value in held["compound_id"]])
        if cfg.get("refit", False):
            if method == "logistic_c0.05":
                model = LogisticRegression(
                    C=0.05, max_iter=4000, solver="liblinear", random_state=20260913 + fold
                )
            elif method == "rbf_svm_c1":
                model = SVC(C=1.0, kernel="rbf")
            elif method == "complement_nb_a1":
                model = ComplementNB(alpha=1.0)
            else:
                raise ValueError(f"Unknown fixed model: {method}")
            valid = train_folds != fold
            model.fit(features[train_indices[valid]], train_y[valid])
        else:
            model = joblib.load(inp / f"model-{method}-fold{fold}.joblib")
        x = features[indices]
        score = (
            model.predict_proba(x)[:, 1]
            if hasattr(model, "predict_proba")
            else model.decision_function(x)
        )
        if not np.isfinite(score).all():
            raise ValueError(f"Nonfinite Morgan control score: {method} {fold}")
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
