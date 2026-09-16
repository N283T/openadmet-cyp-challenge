"""Paired CYP mechanism-feature additions under the retained regression parent CV."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from scipy.stats import spearmanr

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(index: int) -> pl.DataFrame:
    return pl.read_parquet(next(inp.glob(f"{index:03d}-*")))


def st_rae(y: np.ndarray, prediction: np.ndarray, lower: np.ndarray, upper: np.ndarray) -> float:
    def interval_error(values: np.ndarray) -> np.ndarray:
        return np.maximum(values - upper, 0) + np.maximum(lower - values, 0)

    denominator = interval_error(np.full_like(y, y.mean())).sum()
    if denominator <= 0:
        raise ValueError("ST-RAE denominator is not positive")
    return float(interval_error(prediction).sum() / denominator)


def score(frame: pl.DataFrame) -> dict[str, float]:
    y = frame["y_true"].to_numpy()
    prediction = frame["prediction"].to_numpy()
    lower = frame["conf_low"].to_numpy()
    upper = frame["conf_high"].to_numpy()
    return {
        "ST-RAE": st_rae(y, prediction, lower, upper),
        "MAE": float(np.abs(y - prediction).mean()),
        "Spearman": float(spearmanr(y, prediction).statistic),
    }


embedding, labels, assignments, nested = [read(index) for index in range(4)]
feature_tables = {
    "smartcyp": read(4),
    "attnsom_properties": read(5),
    "attnsom_emb": read(6),
    "cypmap_properties": read(7),
    "cypmap_emb": read(8),
}
previous_control = read(9)
embedding_columns = sorted(column for column in embedding.columns if column.startswith("emb_"))
if len(embedding_columns) != 600:
    raise ValueError(f"Expected 600 base embedding columns, found {len(embedding_columns)}")

feature_columns: dict[str, list[str]] = {}
for name, table in feature_tables.items():
    if table["compound_id"].n_unique() != table.height:
        raise ValueError(f"Duplicate compound IDs in {name}")
    feature_columns[name] = [column for column in table.columns if column != "compound_id"]

rows: list[pl.DataFrame] = []
fits: list[dict] = []
memberships: list[pl.DataFrame] = []
endpoints = sorted(
    labels.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))[
        "endpoint"
    ].unique()
)
for endpoint in endpoints:
    data = (
        labels.filter(pl.col("endpoint") == endpoint)
        .join(assignments, on="compound_id", validate="1:1")
        .join(embedding, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    for name, table in feature_tables.items():
        data = data.join(table, on="compound_id", how="left", validate="1:1", suffix=f"_{name}")
    y = data["value"].to_numpy()
    for fold in range(cfg["n_folds"]):
        validation = (data["fold"].to_numpy() == fold) & ~data["is_anchor"].to_numpy()
        auxiliary = data.select("compound_id").join(
            nested.filter(pl.col("outer_fold") == fold),
            on="compound_id",
            how="left",
            validate="1:1",
            maintain_order="left",
        )
        base = np.column_stack(
            [
                data.select(embedding_columns).to_numpy(),
                auxiliary.select([f"aux_{index}" for index in range(4)]).to_numpy(),
            ]
        ).astype(np.float32)
        if not np.isfinite(base).all():
            raise ValueError(f"Non-finite base features for {endpoint} fold {fold}")
        memberships.append(
            data.select("compound_id", "is_anchor", "is_neighbor").with_columns(
                pl.lit(endpoint).alias("endpoint"),
                pl.lit(fold).alias("outer_fold"),
                pl.Series("validation", validation),
            )
        )
        for method, additions in cfg["arms"].items():
            matrices = [base]
            for block in additions:
                matrix = data.select(feature_columns[block]).to_numpy().astype(np.float32)
                matrix[~np.isfinite(matrix)] = np.nan
                matrices.append(matrix)
            features = np.column_stack(matrices)
            model = lgb.LGBMRegressor(**cfg["lightgbm"], random_state=42 + fold)
            model.fit(features[~validation], y[~validation])
            prediction = model.predict(features[validation]).astype(float)
            if not np.isfinite(prediction).all():
                raise ValueError(f"Non-finite predictions for {endpoint} {fold} {method}")
            model.booster_.save_model(str(out / f"{method}-{endpoint}-{fold}.lgbm.txt"))
            rows.append(
                data.filter(pl.Series(validation))
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
                    pl.lit(method).alias("method"),
                    pl.lit(cfg["split_run"]).alias("split_id"),
                    pl.Series("prediction", prediction),
                )
            )
            fits.append(
                {
                    "method": method,
                    "endpoint": endpoint,
                    "fold": fold,
                    "train": int((~validation).sum()),
                    "validation": int(validation.sum()),
                    "dimensions": features.shape[1],
                }
            )
            print(json.dumps(fits[-1]), flush=True)

result = pl.concat(rows)
keys = ["compound_id", "endpoint", "fold"]
control = result.filter(pl.col("method") == "control").sort(keys)
previous_control = previous_control.sort(keys)
if control.height != previous_control.height or control.height != 6384:
    raise ValueError("Control coverage does not match the retained parent readout")
for column in keys + ["y_true", "conf_low", "conf_high"]:
    if not control[column].equals(previous_control[column]):
        raise ValueError(f"Control key/label mismatch in {column}")
max_control_difference = float(
    np.max(np.abs(control["prediction"].to_numpy() - previous_control["prediction"].to_numpy()))
)
if max_control_difference != 0:
    raise ValueError(f"Control prediction reproduction failed: {max_control_difference}")

metric_rows = []
for (endpoint, method), frame in result.group_by("endpoint", "method"):
    metric_rows.append({"endpoint": endpoint, "method": method, "n": frame.height, **score(frame)})
metrics = pl.DataFrame(metric_rows).sort("endpoint", "ST-RAE", "method")
macro = metrics.group_by("method").agg(
    pl.col("ST-RAE").mean().alias("macro_ST-RAE"),
    pl.col("MAE").mean().alias("macro_MAE"),
    pl.col("Spearman").mean().alias("macro_Spearman"),
).sort("macro_ST-RAE")

result.write_parquet(out / "mechanism-oof.parquet")
pl.concat(memberships).write_parquet(out / "memberships.parquet")
metrics.write_csv(out / "metrics.csv")
macro.write_csv(out / "macro.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "embedding_columns": embedding_columns,
            "feature_columns": feature_columns,
            "fits": fits,
            "control_bitwise_equal": True,
            "control_max_abs_difference": max_control_difference,
        },
        indent=2,
    )
    + "\n"
)
print(metrics, flush=True)
print(macro, flush=True)
