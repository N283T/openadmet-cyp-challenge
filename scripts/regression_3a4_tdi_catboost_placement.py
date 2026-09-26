"""Test label-free held-fold placement of a new 3A4 CatBoost peer."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
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
    keys = ["compound_id", "endpoint", "fold"]
    base = (
        pl.read_parquet(inp / "scored-bank-oof.parquet")
        .filter((pl.col("endpoint") == cfg["endpoint"]) & (pl.col("method") == cfg["base_method"]))
        .sort(keys)
    )
    new = pl.read_parquet(inp / "catboost-oof.parquet").sort(keys)
    if base.height != new.height or base.height != 2300:
        raise ValueError("Incomplete aligned 3A4 OOF")
    for column in keys + ["y_true", "conf_low", "conf_high"]:
        if not base[column].equals(new[column]):
            raise ValueError(f"Mismatch: {column}")
    old = base["prediction"].to_numpy()
    raw = new["prediction"].to_numpy()
    folds = base["fold"].to_numpy()
    placed = np.empty_like(raw)
    moments = []
    for fold in range(5):
        mask = folds == fold
        placed[mask] = (raw[mask] - raw[mask].mean()) / raw[mask].std() * old[mask].std() + old[
            mask
        ].mean()
        moments.append(
            {
                "fold": fold,
                "n": int(mask.sum()),
                "raw_mean": float(raw[mask].mean()),
                "raw_sd": float(raw[mask].std()),
                "target_mean": float(old[mask].mean()),
                "target_sd": float(old[mask].std()),
            }
        )
    metrics = [{"method": "scored_bank", "weight": 0.0, "n": base.height, **score(base, old)}]
    for weight in cfg["weights"]:
        for kind, values in [("raw", raw), ("placed", placed)]:
            prediction = (1 - weight) * old + weight * values
            metrics.append(
                {"method": kind, "weight": weight, "n": base.height, **score(base, prediction)}
            )
    pl.DataFrame(metrics).sort("ST_RAE").write_csv(out / "metrics.csv")
    pl.DataFrame(moments).write_csv(out / "moments.csv")
    new.with_columns(pl.Series("placed_prediction", placed)).write_parquet(
        out / "placed-oof.parquet"
    )
    (out / "summary.json").write_text(
        json.dumps(
            {"config": cfg, "label_free_placement": True, "ranks_preserved_within_each_fold": True},
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
