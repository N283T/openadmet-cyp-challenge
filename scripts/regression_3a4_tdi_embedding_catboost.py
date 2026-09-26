"""Train CatBoost readouts on fold-specific TDI continuous embeddings for 3A4."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from catboost import CatBoostRegressor
from scipy.stats import spearmanr


def score(frame: pl.DataFrame, pred: np.ndarray) -> dict[str, float]:
    y = frame["y_true"].to_numpy()
    lo = frame["conf_low"].to_numpy()
    hi = frame["conf_high"].to_numpy()
    center = y.mean()
    denominator = (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    return {
        "ST_RAE": float((np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum() / denominator),
        "MAE": float(np.abs(y - pred).mean()),
        "Spearman": float(spearmanr(y, pred).statistic),
    }


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    endpoint = cfg["endpoint"]
    labels = pl.read_parquet(inp / "measurements.parquet").filter(pl.col("endpoint") == endpoint)
    assignments = pl.read_parquet(inp / "parent-assignments.parquet")
    nested = pl.read_parquet(inp / "nested-log2fc.parquet")
    base = labels.join(assignments, on="compound_id", validate="1:1").sort("compound_id")
    current = (
        pl.read_parquet(inp / "scored-bank-oof.parquet")
        .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == "caruana"))
        .sort("compound_id")
    )
    arm = (
        pl.read_parquet(inp / "tdi-arm-oof.parquet")
        .filter(pl.col("endpoint") == endpoint)
        .sort("compound_id")
    )
    keys = ["compound_id", "endpoint", "fold"]
    rows = []
    fits = []
    for fold in range(5):
        embeddings = pl.read_parquet(inp / f"embedding-fold-{fold}.parquet")
        cols = sorted(column for column in embeddings.columns if column.startswith("emb_"))
        if len(cols) != 600 or embeddings["fold"].unique().to_list() != [fold]:
            raise ValueError(f"Invalid TDI embeddings for fold {fold}")
        if embeddings["split_id"].unique().to_list() != [cfg["split_run"]]:
            raise ValueError("TDI embedding split mismatch")
        features = embeddings.select("compound_id", *cols).join(
            nested.filter(pl.col("outer_fold") == fold).select(
                "compound_id", *[f"aux_{i}" for i in range(4)]
            ),
            on="compound_id",
            validate="1:1",
        )
        data = base.join(features, on="compound_id", validate="1:1").sort("compound_id")
        if data.height != base.height:
            raise ValueError("Incomplete label/feature join")
        columns = cols + [f"aux_{i}" for i in range(4)]
        x = data.select(columns).to_numpy().astype(np.float32)
        if not np.isfinite(x).all():
            raise ValueError("Nonfinite TDI features")
        valid = (data["fold"].to_numpy() == fold) & ~data["is_anchor"].to_numpy()
        model = CatBoostRegressor(**cfg["catboost"], random_seed=42 + fold)
        model.fit(x[~valid], data["value"].to_numpy()[~valid])
        prediction = model.predict(x[valid]).astype(float)
        if not np.isfinite(prediction).all():
            raise ValueError("Nonfinite CatBoost prediction")
        model.save_model(str(out / f"catboost-fold-{fold}.cbm"))
        rows.append(
            data.filter(pl.Series(valid))
            .select(
                "compound_id",
                "endpoint",
                "fold",
                "is_neighbor",
                "conf_low",
                "conf_high",
                pl.col("value").alias("y_true"),
            )
            .with_columns(
                pl.lit("tdi_embedding_catboost").alias("method"),
                pl.Series("prediction", prediction),
            )
        )
        fits.append(
            {
                "fold": fold,
                "train": int((~valid).sum()),
                "valid": int(valid.sum()),
                "features": x.shape[1],
                "seed": 42 + fold,
            }
        )
        print(json.dumps(fits[-1]), flush=True)
    oof = pl.concat(rows).sort(keys)
    if oof.height != 2300 or oof.select(keys).n_unique() != 2300:
        raise ValueError("Incomplete CatBoost OOF")
    for reference in [current, arm]:
        if reference.height != oof.height:
            raise ValueError("Reference OOF coverage mismatch")
        for column in keys + ["y_true", "conf_low", "conf_high"]:
            if not reference[column].equals(oof[column]):
                raise ValueError(f"Reference mismatch: {column}")
    oof.write_parquet(out / "catboost-oof.parquet")
    raw = oof["prediction"].to_numpy()
    anchor = current["prediction"].to_numpy()
    lgbm = arm["prediction"].to_numpy()
    metrics = []
    for name, prediction in [
        ("scored_bank", anchor),
        ("tdi_arm_lgbm", lgbm),
        ("tdi_embedding_catboost", raw),
    ]:
        metrics.append(
            {
                "method": name,
                "blend_weight": 1.0 if name != "scored_bank" else 0.0,
                "n": oof.height,
                **score(oof, prediction),
            }
        )
    for weight in cfg["blend_weights"]:
        if weight == 0:
            continue
        prediction = (1 - weight) * anchor + weight * raw
        metrics.append(
            {
                "method": "scored_plus_catboost",
                "blend_weight": weight,
                "n": oof.height,
                **score(oof, prediction),
            }
        )
    pl.DataFrame(metrics).sort("ST_RAE").write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "fits": fits,
                "classification_representation_fold_safe": True,
                "same_fold_3a4_labels_excluded_from_fit": True,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
