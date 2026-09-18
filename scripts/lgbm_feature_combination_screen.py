"""Screen heterogeneous feature-block combinations with fixed LightGBM readouts."""

import importlib.util
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


def read(name: str) -> pl.DataFrame:
    path = inp / name
    if not path.is_file():
        raise FileNotFoundError(path)
    return pl.read_parquet(path)


spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)

assignments = read("parent-assignments.parquet")
measurements = read("measurements.parquet").filter(
    pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
)

blocks = {}
block_columns = {}
fold_blocks = set()
for block_name, block_cfg in cfg["blocks"].items():
    frame = read(block_cfg["file"])
    source_columns = sorted(
        column
        for column in frame.columns
        if any(column.startswith(prefix) for prefix in block_cfg["prefixes"])
    )
    if not source_columns:
        raise ValueError(f"No columns for block {block_name}")
    renamed = {column: f"{block_name}::{column}" for column in source_columns}
    keep = ["compound_id", *source_columns]
    if "outer_fold_column" in block_cfg:
        keep.insert(1, block_cfg["outer_fold_column"])
        fold_blocks.add(block_name)
    blocks[block_name] = frame.select(keep).rename(renamed)
    block_columns[block_name] = list(renamed.values())


def make_model(fold: int) -> lgb.LGBMRegressor:
    params = dict(cfg["lightgbm"])
    params.update(n_jobs=cfg["num_threads"], random_state=cfg["seed"] + fold)
    return lgb.LGBMRegressor(**params)


def score(rows: pl.DataFrame, prediction: np.ndarray) -> float:
    return float(
        scorer.rae_soft_threshold_absolute_error(
            rows["y_true"].to_numpy(),
            prediction,
            y_true_lower=rows["conf_low"].to_numpy(),
            y_true_upper=rows["conf_high"].to_numpy(),
        )
    )


def matrix_for(candidate: dict, table: pl.DataFrame, fold: int):
    matrices = []
    names = []
    labels = []
    ids = table.select("compound_id")
    for block_name in candidate["blocks"]:
        frame = blocks[block_name]
        if block_name in fold_blocks:
            fold_column = cfg["blocks"][block_name]["outer_fold_column"]
            frame = frame.filter(pl.col(fold_column) == fold).drop(fold_column)
        aligned = ids.join(
            frame, on="compound_id", how="left", validate="1:1", maintain_order="left"
        )
        columns = block_columns[block_name]
        if aligned.select(columns).null_count().sum_horizontal()[0] and block_name != "mordred":
            raise ValueError(f"Unexpected nulls in {block_name}, fold {fold}")
        values = aligned.select(columns).to_numpy().astype(np.float32)
        values[~np.isfinite(values)] = np.nan
        matrices.append(values)
        names.extend(columns)
        labels.extend([block_name] * len(columns))
    return np.column_stack(matrices), np.asarray(names), np.asarray(labels)


oof_rows = []
metric_rows = []
fit_rows = []
importance_rows = []
selection_rows = []
for candidate in cfg["candidates"]:
    endpoint_rows = []
    for endpoint in sorted(measurements["endpoint"].unique()):
        table = (
            measurements.filter(pl.col("endpoint") == endpoint)
            .join(assignments, on="compound_id", validate="1:1")
            .sort("compound_id")
            .rename({"value": "y_true"})
        )
        y = table["y_true"].to_numpy()
        folds = table["fold"].to_numpy()
        anchors = table["is_anchor"].to_numpy()
        for fold in range(5):
            validation = (folds == fold) & ~anchors
            training = ~validation
            x, feature_names, feature_blocks = matrix_for(candidate, table, fold)
            selected = np.arange(x.shape[1])
            if "importance_top_k" in candidate:
                selector = make_model(fold)
                selector.fit(x[training], y[training])
                gain = selector.booster_.feature_importance(importance_type="gain")
                order = np.lexsort((feature_names, -gain))
                selected = order[: min(candidate["importance_top_k"], len(order))]
                for rank, index in enumerate(selected, start=1):
                    selection_rows.append(
                        {
                            "method": candidate["name"],
                            "endpoint": endpoint,
                            "fold": fold,
                            "rank": rank,
                            "feature": feature_names[index],
                            "block": feature_blocks[index],
                            "selector_gain": float(gain[index]),
                        }
                    )
            model = make_model(fold)
            model.fit(x[training][:, selected], y[training])
            prediction = np.asarray(
                model.predict(x[validation][:, selected]), dtype=float
            )
            gain = model.booster_.feature_importance(importance_type="gain")
            split = model.booster_.feature_importance(importance_type="split")
            for index, feature_index in enumerate(selected):
                importance_rows.append(
                    {
                        "method": candidate["name"],
                        "endpoint": endpoint,
                        "fold": fold,
                        "feature": feature_names[feature_index],
                        "block": feature_blocks[feature_index],
                        "gain": float(gain[index]),
                        "split": int(split[index]),
                    }
                )
            rows = table.filter(pl.Series(validation)).select(
                "compound_id",
                "endpoint",
                "fold",
                "is_neighbor",
                "conf_low",
                "conf_high",
                "y_true",
            ).with_columns(
                pl.lit(candidate["name"]).alias("method"),
                pl.lit(cfg["split_run"]).alias("split_id"),
                pl.Series("prediction", prediction),
            )
            endpoint_rows.append(rows)
            fit = {
                "method": candidate["name"],
                "endpoint": endpoint,
                "fold": fold,
                "input_dimensions": x.shape[1],
                "selected_dimensions": len(selected),
                "train": int(training.sum()),
                "validation": int(validation.sum()),
                "fold_ST_RAE": score(rows, prediction),
            }
            fit_rows.append(fit)
            print(json.dumps(fit), flush=True)
        endpoint_oof = pl.concat(endpoint_rows[-5:])
        prediction = endpoint_oof["prediction"].to_numpy()
        y_true = endpoint_oof["y_true"].to_numpy()
        metric_rows.append(
            {
                "method": candidate["name"],
                "endpoint": endpoint,
                "ST_RAE": score(endpoint_oof, prediction),
                "MAE": float(np.abs(y_true - prediction).mean()),
                "Spearman": float(spearmanr(y_true, prediction).statistic),
                "n": endpoint_oof.height,
            }
        )
    candidate_oof = pl.concat(endpoint_rows).sort("endpoint", "fold", "compound_id")
    candidate_oof.write_parquet(out / f"{candidate['name']}-oof.parquet")
    oof_rows.append(candidate_oof)

metrics = pl.DataFrame(metric_rows).sort("endpoint", "ST_RAE")
macro = metrics.group_by("method").agg(
    pl.col("ST_RAE", "MAE", "Spearman").mean()
).sort("ST_RAE")
fits = pl.DataFrame(fit_rows)
importance = pl.DataFrame(importance_rows)
block_importance = (
    importance.group_by("method", "endpoint", "fold", "block")
    .agg(pl.col("gain", "split").sum())
    .sort("method", "endpoint", "fold", "gain", descending=[False, False, False, True])
)
metrics.write_csv(out / "metrics.csv")
macro.write_csv(out / "macro.csv")
fits.write_csv(out / "fits.csv")
importance.write_parquet(out / "feature-importance.parquet")
block_importance.write_csv(out / "block-importance.csv")
if selection_rows:
    pl.DataFrame(selection_rows).write_parquet(out / "selected-features.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "primary_selection_scope": "endpoint",
            "macro_use": "descriptive_only",
            "feature_selection": (
                "Outer-training-only LightGBM gain ranking, then refit on the same "
                "outer-training rows."
            ),
            "endpoint_metrics": metrics.to_dicts(),
            "macro": macro.to_dicts(),
        },
        indent=2,
    )
    + "\n"
)
print(macro, flush=True)
print(metrics, flush=True)
