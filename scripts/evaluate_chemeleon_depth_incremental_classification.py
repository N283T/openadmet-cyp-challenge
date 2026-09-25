"""Labelled-OOF ranking diagnostic for depth additions to submitted TDI models."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    ranked = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        ranked[mask] = rankdata(values[mask], method="average") / (mask.sum() + 1)
    return ranked


def main():
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    keys = ["compound_id", "endpoint", "fold"]
    d6 = pl.read_parquet(inp / cfg["base_d6"]).select(
        *keys, "y_true", pl.col("prediction").alias("base")
    )
    a4 = (
        pl.read_parquet(inp / cfg["base_a4"])
        .filter(pl.col("model") == cfg["a4_model"])
        .with_columns(pl.lit("CYP3A4_is_TDI").alias("endpoint"))
        .select(*keys, "y_true", pl.col("prediction").alias("base"))
    )
    base = pl.concat([d6, a4]).sort(keys)
    depth = pl.read_parquet(inp / cfg["depth_oof"]).filter(pl.col("track") == "classification")
    assert depth.height == 3 * base.height
    for level in (1, 3, 6):
        single = depth.filter(pl.col("depth") == level).select(
            *keys, pl.col("y_true").alias("depth_y"), pl.col("prediction").alias(f"d{level}")
        )
        base = base.join(single, on=keys, validate="1:1")
        assert base.height == 5081
        assert np.array_equal(base["y_true"].to_numpy(), base["depth_y"].to_numpy())
        base = base.drop("depth_y")
    rows = []
    for endpoint, part in base.group_by("endpoint", maintain_order=True):
        endpoint = endpoint[0]
        y = part["y_true"].to_numpy().astype(int)
        folds = part["fold"].to_numpy()
        base_score = fold_rank(part["base"].to_numpy(), folds)
        for levels in cfg["depth_sets"]:
            depth_score = fold_rank(
                np.mean([part[f"d{level}"].to_numpy() for level in levels], axis=0), folds
            )
            corr = float(spearmanr(base_score, depth_score).statistic)
            for weight in cfg["weights"]:
                score = (1 - weight) * base_score + weight * depth_score
                rows.append(
                    {
                        "endpoint": endpoint,
                        "depths": "+".join(map(str, levels)),
                        "weight": weight,
                        "n": len(y),
                        "base_depth_spearman": corr,
                        "AUROC": float(roc_auc_score(y, score)),
                        "AP": float(average_precision_score(y, score)),
                    }
                )
    metrics = pl.DataFrame(rows).sort("endpoint", "depths", "weight")
    metrics.write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "description": (
                    "Labelled-only OOF rank screen; deployment needs complete-pool "
                    "D6/A4 refit and call-rate evaluation"
                ),
                "config": cfg,
            },
            indent=2,
        )
        + "\n"
    )
    print(metrics, flush=True)


if __name__ == "__main__":
    main()
