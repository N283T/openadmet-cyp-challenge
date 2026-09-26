"""Audit adding a standalone model or complete pair block to ten-member 1A2."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr


def st_rae(y: np.ndarray, lo: np.ndarray, hi: np.ndarray, pred: np.ndarray) -> float:
    error = (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
    mean = y.mean()
    baseline = (np.maximum(lo - mean, 0) + np.maximum(mean - hi, 0)).sum()
    return float(error / baseline)


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    keys = ["compound_id", "endpoint", "fold"]
    current = pl.read_parquet(inp / "current-oof.parquet").filter(
        (pl.col("endpoint") == cfg["endpoint"]) & (pl.col("scenario") == "current")
    ).select(*keys, pl.col("prediction").alias("current"))
    atomic = pl.read_parquet(inp / "atomic-oof.parquet").filter(
        (pl.col("endpoint") == cfg["endpoint"]) & (pl.col("split_id") == "parent")
        & pl.col("method").is_in(cfg["pair_members"])
    )
    wide = atomic.pivot(
        on="method", index=keys + ["y_true", "conf_low", "conf_high"], values="prediction"
    ).join(current, on=keys, validate="1:1").sort("compound_id")
    if wide.height != 1372 or current.height != 1372:
        raise ValueError("Incomplete current/atomic OOF alignment")
    control = wide["current"].to_numpy()
    jeremy = wide["jeremy_anchor_lgbm"].to_numpy()
    pfn = wide["all_pool_tabpfn"].to_numpy()
    count = cfg["existing_members"]
    predictions = {
        "current": control,
        "add_jeremy_anchor": (count * control + jeremy) / (count + 1),
        "add_pair_block": (count * control + (jeremy + pfn) / 2) / (count + 1),
    }
    rows = []
    outputs = []
    for name, pred in predictions.items():
        outputs.append(wide.select(*keys).with_columns(
            pl.lit(name).alias("scenario"), pl.Series("prediction", pred)
        ))
        for fold in [-1, 0, 1, 2, 3, 4]:
            mask = (
                np.ones(wide.height, dtype=bool)
                if fold == -1 else wide["fold"].to_numpy() == fold
            )
            y = wide["y_true"].to_numpy()[mask]
            lo = wide["conf_low"].to_numpy()[mask]
            hi = wide["conf_high"].to_numpy()[mask]
            p = pred[mask]
            rows.append({
                "scenario": name, "fold": fold, "n": len(y),
                "ST_RAE": st_rae(y, lo, hi, p),
                "MAE": float(np.abs(y - p).mean()),
                "Spearman": float(spearmanr(y, p).statistic),
            })
    pl.DataFrame(rows).sort("scenario", "fold").write_csv(out / "metrics.csv")
    pl.concat(outputs).write_parquet(out / "oof.parquet")
    (out / "summary.json").write_text(json.dumps({
        "config": cfg, "rows": wide.height,
        "qualification": "Parent-only incremental check; not a scaffold or test result",
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
