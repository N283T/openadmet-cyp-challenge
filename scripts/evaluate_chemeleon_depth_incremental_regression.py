"""Fixed-weight OOF addition to the currently submitted regression components."""

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
    keys = ["compound_id", "endpoint", "fold"]
    endpoints = {
        "CYP1A2": "base_1a2",
        "CYP2C9": "base_2c9",
        "CYP3A4": "base_3a4",
    }
    frames = []
    for endpoint, source in endpoints.items():
        frame = pl.read_parquet(inp / cfg[source]).filter(
            (pl.col("endpoint") == f"{endpoint}_pIC50_direct_inhibition")
            & (pl.col("method") == cfg["methods"][source])
        )
        frames.append(
            frame.select(
                *keys, "y_true", "conf_low", "conf_high", pl.col("prediction").alias("base")
            )
        )
    tabpfn = pl.read_parquet(inp / cfg["base_d6_tabpfn"]).filter(
        pl.col("endpoint") == "CYP2D6_pIC50_direct_inhibition"
    )
    lgbm = pl.read_parquet(inp / cfg["base_d6_lgbm"]).filter(
        pl.col("endpoint") == "CYP2D6_pIC50_direct_inhibition"
    )
    d6 = tabpfn.select(
        *keys, "y_true", "conf_low", "conf_high", pl.col("prediction").alias("tabpfn")
    )
    d6 = d6.join(lgbm.select(*keys, pl.col("prediction").alias("lgbm")), on=keys, validate="1:1")
    d6 = d6.with_columns(((pl.col("tabpfn") + pl.col("lgbm")) / 2).alias("base"))
    frames.append(d6.select(*keys, "y_true", "conf_low", "conf_high", "base"))
    base = pl.concat(frames).sort(keys)
    assert base.height == 6384 and base.select(keys).unique().height == base.height
    depth = pl.read_parquet(inp / cfg["depth_oof"]).filter(pl.col("track") == "regression")
    assert depth.height == 3 * base.height
    for level in (1, 3, 6):
        single = depth.filter(pl.col("depth") == level).select(
            *keys, pl.col("y_true").alias("depth_y"), pl.col("prediction").alias(f"d{level}")
        )
        base = base.join(single, on=keys, validate="1:1")
        assert base.height == 6384
        assert np.allclose(base["y_true"].to_numpy(), base["depth_y"].to_numpy())
        base = base.drop("depth_y")
    rows = []
    for endpoint, part in base.group_by("endpoint", maintain_order=True):
        endpoint = endpoint[0]
        y = part["y_true"].to_numpy()
        low = part["conf_low"].to_numpy()
        high = part["conf_high"].to_numpy()
        base_prediction = part["base"].to_numpy()
        folds = part["fold"].to_numpy()
        for levels in cfg["depth_sets"]:
            depth_prediction = np.mean([part[f"d{level}"].to_numpy() for level in levels], axis=0)
            corr = float(spearmanr(base_prediction, depth_prediction).statistic)
            for weight in cfg["weights"]:
                prediction = (1 - weight) * base_prediction + weight * depth_prediction
                for fold in [-1, 0, 1, 2, 3, 4]:
                    mask = np.ones(len(y), dtype=bool) if fold == -1 else folds == fold
                    rows.append(
                        {
                            "endpoint": endpoint,
                            "depths": "+".join(map(str, levels)),
                            "weight": weight,
                            "fold": fold,
                            "n": int(mask.sum()),
                            "base_depth_spearman": corr,
                            "ST_RAE": float(
                                scorer.rae_soft_threshold_absolute_error(
                                    y[mask],
                                    prediction[mask],
                                    y_true_lower=low[mask],
                                    y_true_upper=high[mask],
                                )
                            ),
                            "MAE": float(np.mean(np.abs(y[mask] - prediction[mask]))),
                        }
                    )
    metrics = pl.DataFrame(rows).sort("endpoint", "depths", "weight", "fold")
    metrics.write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "description": (
                    "Fixed, untuned incremental OOF blends with current submitted "
                    "regression component recipes"
                ),
                "n": base.height,
                "config": cfg,
            },
            indent=2,
        )
        + "\n"
    )
    print(
        metrics.filter(pl.col("fold") == -1).select("endpoint", "depths", "weight", "ST_RAE"),
        flush=True,
    )


if __name__ == "__main__":
    main()
