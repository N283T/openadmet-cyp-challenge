"""Compare complete D6 ensemble members and internal equal-depth blocks."""

import importlib.util
import json
import os
from pathlib import Path

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
    endpoint = "CYP2D6_pIC50_direct_inhibition"
    keys = ["compound_id", "endpoint", "fold"]
    base = None
    for name, source in cfg["sources"].items():
        frame = pl.read_parquet(inp / source["file"]).filter(pl.col("endpoint") == endpoint)
        if "depth" in source:
            frame = frame.filter(pl.col("depth") == source["depth"])
        assert frame.height == 1457 and frame.select(keys).unique().height == 1457, name
        if base is None:
            base = frame.select(*keys, "y_true", "conf_low", "conf_high").sort(keys)
        aligned = frame.select(*keys, pl.col("prediction").alias(name))
        base = base.join(aligned, on=keys, validate="1:1").sort(keys)
        assert base.height == 1457
        assert np.allclose(base["y_true"].to_numpy(), frame.sort(keys)["y_true"].to_numpy())
    for name, members in cfg["blocks"].items():
        base = base.with_columns(pl.mean_horizontal(*members).alias(name))
    y = base["y_true"].to_numpy()
    low = base["conf_low"].to_numpy()
    high = base["conf_high"].to_numpy()
    folds = base["fold"].to_numpy()
    control = (
        base["current_tabpfn"].to_numpy() + base["current_ft_lgbm"].to_numpy()
    ) / 2
    rows = []
    prediction_frames = []
    recipes = {name: [name] for name in cfg["sources"] | cfg["blocks"]}
    recipes.update(cfg["scenarios"])
    for name, members in recipes.items():
        pred = np.mean([base[member].to_numpy() for member in members], axis=0)
        prediction_frames.append(base.select(*keys).with_columns(
            pl.lit(name).alias("scenario"), pl.Series("prediction", pred)
        ))
        for fold in [-1, 0, 1, 2, 3, 4]:
            mask = np.ones(len(y), dtype=bool) if fold == -1 else folds == fold
            rows.append({
                "scenario": name,
                "kind": "ensemble" if name in cfg["scenarios"] else "single_or_block",
                "members": "+".join(members),
                "fold": fold,
                "n": int(mask.sum()),
                "ST_RAE": float(scorer.rae_soft_threshold_absolute_error(
                    y[mask], pred[mask], y_true_lower=low[mask], y_true_upper=high[mask]
                )),
                "MAE": float(np.abs(y[mask] - pred[mask]).mean()),
                "Spearman_truth": float(spearmanr(y[mask], pred[mask]).statistic),
                "Spearman_control": float(spearmanr(control[mask], pred[mask]).statistic),
            })
    pl.DataFrame(rows).sort("kind", "scenario", "fold").write_csv(out / "metrics.csv")
    pl.concat(prediction_frames).write_parquet(out / "oof.parquet")
    (out / "summary.json").write_text(json.dumps({
        "config": cfg,
            "rule": (
                "Each depth block averages its complete member predictions; "
                "the outer ensemble averages its complete members equally. "
                "No fractional mixing sweep."
            ),
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
