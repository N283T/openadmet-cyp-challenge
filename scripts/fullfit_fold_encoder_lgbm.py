"""Fit one full-label LightGBM per fold-specific encoder and average test predictions."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl


def main():
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    labels = pl.read_parquet(inp / cfg["measurements"]).filter(
        pl.col("endpoint") == cfg["endpoint"]
    )
    assignments = pl.read_parquet(inp / cfg["assignments"])
    aliases = (
        pl.read_parquet(inp / cfg["aliases"])
        .filter(pl.col("role") == "test")
        .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
        .unique()
    )
    reference = pl.read_csv(inp / cfg["reference"])
    test = reference.select("SMILES", "Molecule_Name").join(
        aliases, on="Molecule_Name", validate="1:1", maintain_order="left"
    )
    auxiliary = pl.read_parquet(inp / cfg["auxiliary"])
    auxcols = [f"aux_{i}" for i in range(4)]
    assert labels.height == labels["compound_id"].n_unique()
    assert set(labels["compound_id"]).issubset(set(assignments["compound_id"]))
    assert test.height == test["compound_id"].n_unique() == 750
    assert auxiliary["compound_id"].n_unique() == auxiliary.height == 6896

    fold_rows = []
    fits = []
    feature_files = cfg["feature_files"]
    assert sorted(map(int, feature_files)) == list(range(5))
    for fold in range(5):
        features = pl.read_parquet(inp / feature_files[str(fold)])
        embcols = sorted(c for c in features.columns if c.startswith(cfg["feature_prefix"]))
        assert len(embcols) == cfg["embedding_dimensions"]
        assert features["compound_id"].n_unique() == features.height == 6896
        features = features.select("compound_id", *embcols).join(
            auxiliary.select("compound_id", *auxcols),
            on="compound_id",
            validate="1:1",
        )
        cols = embcols + auxcols
        train = labels.join(features, on="compound_id", validate="1:1").sort("compound_id")
        blind = test.join(features, on="compound_id", validate="1:1", maintain_order="left")
        x = train.select(cols).to_numpy().astype(np.float32)
        xt = blind.select(cols).to_numpy().astype(np.float32)
        y = train["value"].to_numpy()
        assert np.isfinite(x).all() and np.isfinite(xt).all() and np.isfinite(y).all()
        model = lgb.LGBMRegressor(
            objective="regression_l1",
            n_estimators=300,
            num_leaves=15,
            min_child_samples=20,
            learning_rate=0.05,
            reg_lambda=1.0,
            n_jobs=4,
            random_state=42 + fold,
            verbosity=-1,
            deterministic=True,
            force_col_wise=True,
        )
        model.fit(x, y)
        prediction = np.asarray(model.predict(xt), dtype=float)
        assert np.isfinite(prediction).all()
        model.booster_.save_model(str(out / f"fold{fold}.lgbm.txt"))
        fold_rows.append(
            blind.select("compound_id", "Molecule_Name", "SMILES").with_columns(
                pl.lit(cfg["endpoint"]).alias("endpoint"),
                pl.lit(fold).alias("fold"),
                pl.lit(cfg["method"]).alias("method"),
                pl.Series("prediction", prediction),
            )
        )
        fits.append(
            {
                "fold_encoder": fold,
                "train_n": train.height,
                "test_n": blind.height,
                "dimensions": len(cols),
                "prediction_mean": float(prediction.mean()),
                "prediction_std": float(prediction.std()),
            }
        )
        print(json.dumps(fits[-1]), flush=True)

    fold_test = pl.concat(fold_rows)
    fold_test.write_parquet(out / "fold-test-predictions.parquet")
    test_predictions = (
        fold_test.group_by("compound_id", "Molecule_Name", "SMILES", "endpoint", "method")
        .agg(pl.col("prediction").mean(), pl.len().alias("model_count"))
        .with_columns(pl.lit(-1).alias("fold"), pl.lit("full_train").alias("split_id"))
    )
    assert test_predictions.height == 750 and (test_predictions["model_count"] == 5).all()
    test_predictions.write_parquet(out / "test-predictions.parquet")
    labels.select("compound_id", "endpoint", "value").with_columns(
        pl.lit("full_training").alias("role")
    ).write_parquet(out / "membership.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "fits": fits,
                "deployment": (
                    "Equal mean of five full-label LGBM readouts over fold-specific encoders"
                ),
                "test_labels_used": False,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
