"""Refit the current FT600 LightGBM member on each Jeremy message-passing depth."""

import importlib.util
import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from scipy.stats import spearmanr


def load_scorer(path):
    spec = importlib.util.spec_from_file_location("official_scoring", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    scorer = load_scorer(inp / cfg["scorer"])
    assignments = pl.read_parquet(inp / cfg["assignments"]).filter(pl.col("strategy") == "series")
    labels = pl.read_parquet(inp / cfg["measurements"]).filter(
        pl.col("endpoint") == cfg["endpoint"]
    )
    nested = pl.read_parquet(inp / cfg["nested"])
    original = pl.read_parquet(inp / cfg["original_oof"]).filter(
        pl.col("endpoint") == cfg["endpoint"]
    )
    columns = [f"emb_{i:03d}" for i in range(cfg["dimensions"])]
    auxiliary = [f"aux_{i}" for i in range(4)]
    predictions = []
    for fold in range(5):
        selected = json.loads((inp / f"fold{fold}-selected.json").read_text())
        outer = nested.filter(pl.col("outer_fold") == fold).select("compound_id", *auxiliary)
        for depth in cfg["depths"]:
            source = (
                f"fold{fold}-original-final.parquet"
                if depth == 4
                else f"fold{fold}-depth{depth}.parquet"
            )
            feature = pl.read_parquet(inp / source)
            assert feature["fold"].unique().to_list() == [fold]
            assert feature["split_id"].unique().to_list() == [cfg["split_id"]]
            table = (
                labels.join(assignments, on="compound_id", validate="1:1")
                .join(feature.select("compound_id", *columns), on="compound_id", validate="1:1")
                .join(outer, on="compound_id", validate="1:1")
                .sort("compound_id")
            )
            validation = (table["fold"].to_numpy() == fold) & ~table["is_anchor"].to_numpy()
            x = table.select(columns + auxiliary).to_numpy().astype(np.float32)
            y = table["value"].to_numpy()
            assert np.isfinite(x).all() and validation.sum() > 0
            model = lgb.LGBMRegressor(
                **selected["params"],
                n_jobs=3,
                random_state=42 + fold,
                verbosity=-1,
                deterministic=True,
                force_col_wise=True,
            )
            model.fit(x[~validation], y[~validation])
            prediction = model.predict(x[validation])
            predictions.append(
                table.filter(pl.Series(validation))
                .select("compound_id", "endpoint", "fold", "conf_low", "conf_high")
                .with_columns(
                    pl.Series("y_true", y[validation]),
                    pl.Series("prediction", prediction),
                    pl.lit(depth).alias("depth"),
                    pl.lit(cfg["split_id"]).alias("split_id"),
                )
            )
            print(f"fold={fold} depth={depth} n={validation.sum()}", flush=True)
    oof = pl.concat(predictions).sort("depth", "compound_id")
    assert oof.height == 3 * 1457
    final = oof.filter(pl.col("depth") == 4).sort("compound_id")
    original = original.sort("compound_id")
    assert final.select("compound_id", "fold", "y_true").equals(
        original.select("compound_id", "fold", "y_true")
    )
    max_difference = float(
        np.max(np.abs(final["prediction"].to_numpy() - original["prediction"].to_numpy()))
    )
    assert max_difference < 1e-5, max_difference
    rows = []
    for (depth,), part in oof.group_by("depth"):
        y = part["y_true"].to_numpy()
        p = part["prediction"].to_numpy()
        rows.append(
            {
                "depth": depth,
                "n": part.height,
                "ST_RAE": float(
                    scorer.rae_soft_threshold_absolute_error(
                        y,
                        p,
                        y_true_lower=part["conf_low"].to_numpy(),
                        y_true_upper=part["conf_high"].to_numpy(),
                    )
                ),
                "MAE": float(np.abs(y - p).mean()),
                "Spearman": float(spearmanr(y, p).statistic),
            }
        )
    oof.write_parquet(out / "oof.parquet")
    pl.DataFrame(rows).sort("depth").write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "final_depth_max_abs_vs_original": max_difference,
                "selected_hyperparameters": (
                    "Reused unchanged from original fold-specific inner tuning"
                ),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
