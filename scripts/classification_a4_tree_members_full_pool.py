"""Extend the two submitted A4 LightGBM members to the complete UMAP pool."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator


def main() -> None:
    inputs = Path(os.environ["CYP_INPUT_DIR"])
    outputs = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    RDLogger.DisableLog("rdApp.*")
    compounds = pl.read_parquet(inputs / cfg["compounds"]).sort("compound_id")
    assignments = (
        pl.read_parquet(inputs / cfg["assignments"])
        .filter(pl.col("strategy") == cfg["split_id"])
        .select("compound_id", "fold", "group")
        .sort("compound_id")
    )
    labels = (
        pl.read_parquet(inputs / cfg["labels"])
        .filter(pl.col("endpoint") == cfg["endpoint"])
        .select("compound_id", pl.col("value").cast(pl.Int64).alias("y_true"))
    )
    pool = assignments.join(labels, on="compound_id", how="left", validate="1:1")
    if pool.height != cfg["expected_pool_rows"]:
        raise ValueError("Incomplete A4 molecular pool")
    ids = compounds.select("compound_id")
    id_to_index = {value: row for row, value in enumerate(ids["compound_id"])}
    emb = pl.read_parquet(inputs / cfg["class_embeddings"])
    class_columns = [column for column in emb.columns if column.startswith("emb_")]
    x_class = (
        ids.join(emb, on="compound_id", validate="1:1")
        .select(class_columns)
        .to_numpy()
        .astype(np.float32)
    )
    mordred = pl.concat([pl.read_parquet(p) for p in sorted(inputs.glob(cfg["mordred_files"]))])
    mordred_columns = [column for column in mordred.columns if column != "compound_id"]
    x_mordred = (
        ids.join(mordred, on="compound_id", validate="1:1")
        .select(mordred_columns)
        .to_numpy()
        .astype(np.float32)
    )
    x_mordred[~np.isfinite(x_mordred)] = np.nan
    x_class_mordred = np.column_stack([x_class, x_mordred])
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    x_morgan = np.stack(
        [
            generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(s))
            for s in compounds["standardized_smiles"]
        ]
    ).astype(np.float32)
    pxr_table = pl.read_parquet(inputs / cfg["pxr_features"])
    pxr_columns = [c for c in pxr_table.columns if c.startswith("pxr_")]
    pxr = (
        ids.join(pxr_table, on="compound_id", validate="1:1")
        .select(pxr_columns)
        .to_numpy()
        .astype(np.float32)
    )
    admetica_table = pl.read_parquet(inputs / cfg["admetica_features"])
    admetica_columns = [c for c in admetica_table.columns if c.endswith("_probability")]
    admetica = (
        ids.join(admetica_table, on="compound_id", validate="1:1")
        .select(admetica_columns)
        .to_numpy()
        .astype(np.float32)
    )
    x_pxr = np.column_stack([x_morgan, pxr, admetica])
    saved = {
        "external_class_mordred": pl.read_parquet(inputs / cfg["lgbm_oof"]),
        "pxr_admetica": pl.read_parquet(inputs / cfg["pxr_oof"]),
    }
    parts = []
    audits = []
    for fold in range(5):
        held = pool.filter(pl.col("fold") == fold)
        held_indices = np.asarray([id_to_index[value] for value in held["compound_id"]])
        labelled = pool.filter(pl.col("y_true").is_not_null())
        train = labelled.filter(pl.col("fold") != fold)
        train_indices = np.asarray([id_to_index[value] for value in train["compound_id"]])
        y_train = train["y_true"].to_numpy().astype(int)
        frozen = lgb.Booster(model_file=str(inputs / f"external-class-mordred-fold{fold}.txt"))
        p_class = frozen.predict(x_class_mordred[held_indices])
        params = {
            "n_estimators": 500,
            "learning_rate": 0.03,
            "num_leaves": 31,
            "min_child_samples": 20,
            "subsample": 0.8,
            "colsample_bytree": 0.6,
            "reg_lambda": 1.0,
            "verbosity": -1,
            "n_jobs": -1,
            "random_state": cfg["pxr_random_seed"] + fold,
        }
        fitted = lgb.LGBMClassifier(**params).fit(x_pxr[train_indices], y_train)
        p_pxr = fitted.predict_proba(x_pxr[held_indices])[:, 1]
        for name, values, method in [
            ("external_class_mordred", p_class, "external_class_mordred_balanced"),
            ("pxr_admetica", p_pxr, "morgan_pxr_admetica_lgbm"),
        ]:
            result = held.with_columns(
                pl.lit(cfg["endpoint"]).alias("endpoint"),
                pl.lit(cfg["split_id"]).alias("split_id"),
                pl.lit(name).alias("model"),
                pl.Series("score", np.asarray(values, dtype=float)),
            )
            expected = saved[name].filter(
                (pl.col("endpoint") == cfg["endpoint"])
                & (pl.col("method") == method)
                & (pl.col("fold") == fold)
            )
            if "strategy" in expected.columns:
                expected = expected.filter(pl.col("strategy") == cfg["split_id"])
            labelled_result = result.filter(pl.col("y_true").is_not_null()).select(
                "compound_id", "score"
            )
            aligned = labelled_result.join(
                expected.select("compound_id", pl.col("prediction").alias("expected")),
                on="compound_id",
                validate="1:1",
            )
            if aligned.height != labelled_result.height:
                raise ValueError(f"Incomplete {name} audit in fold {fold}")
            difference = np.abs(aligned["score"].to_numpy() - aligned["expected"].to_numpy())
            audits.append(
                {
                    "model": name,
                    "fold": fold,
                    "labelled_rows": labelled_result.height,
                    "max_abs_diff": float(difference.max()),
                }
            )
            if difference.max() > 1e-5:
                raise ValueError(f"Replay mismatch in {name} fold {fold}: {difference.max()}")
            parts.append(result)
        print(json.dumps(audits[-2:]), flush=True)
    full = pl.concat(parts).sort("model", "compound_id")
    if full.height != pool.height * 2:
        raise ValueError("Incomplete tree-member OOF")
    full.write_parquet(outputs / "all-pool-oof.parquet")
    (outputs / "audit.json").write_text(json.dumps(audits, indent=2) + "\n")


if __name__ == "__main__":
    main()
