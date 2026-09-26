"""Evaluate fixed equal-component CYP1A2 combinations on matched folds."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr


def st_rae(y: np.ndarray, lo: np.ndarray, hi: np.ndarray, pred: np.ndarray) -> float:
    error = (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
    center = y.mean()
    baseline = (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    return float(error / baseline)


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    keys = ["compound_id", "endpoint", "fold"]
    predictions = []
    for split in ["parent", "scaffold"]:
        aligned = None
        for component in cfg["components"]:
            frame = pl.read_parquet(inp / f"{split}-{component}.parquet").filter(
                pl.col("endpoint") == cfg["endpoint"]
            )
            if frame.height != 1372 or frame.select(keys).is_duplicated().any():
                raise ValueError(f"Unexpected {split}/{component} OOF coverage")
            frame = frame.select(*keys, "y_true", "conf_low", "conf_high",
                                 pl.col("prediction").alias(component))
            if aligned is None:
                aligned = frame
            else:
                aligned = aligned.join(frame, on=keys, validate="1:1", suffix="_new")
                for column in ["y_true", "conf_low", "conf_high"]:
                    if not np.allclose(aligned[column], aligned[f"{column}_new"]):
                        raise ValueError(f"{split}/{component} {column} mismatch")
                    aligned = aligned.drop(f"{column}_new")
        if aligned is None or aligned.height != 1372:
            raise ValueError(f"Incomplete {split} alignment")
        for name, members in cfg["recipes"].items():
            pred = np.mean([aligned[member].to_numpy() for member in members], axis=0)
            predictions.append(aligned.select(*keys, "y_true", "conf_low", "conf_high")
                               .with_columns(pl.lit(split).alias("split_id"),
                                             pl.lit(name).alias("method"),
                                             pl.Series("prediction", pred)))
    oof = pl.concat(predictions).sort("split_id", "method", "compound_id")
    expected = 2 * len(cfg["recipes"]) * 1372
    if oof.height != expected or oof.select("split_id", "method", *keys).is_duplicated().any():
        raise ValueError("Incomplete combination OOF")
    metrics = []
    for (split, method), frame in oof.group_by("split_id", "method"):
        for fold in [-1, 0, 1, 2, 3, 4]:
            part = frame if fold == -1 else frame.filter(pl.col("fold") == fold)
            y = part["y_true"].to_numpy()
            lo = part["conf_low"].to_numpy()
            hi = part["conf_high"].to_numpy()
            pred = part["prediction"].to_numpy()
            metrics.append({
                "split_id": split, "method": method, "fold": fold, "n": part.height,
                "ST_RAE": st_rae(y, lo, hi, pred),
                "MAE": float(np.abs(y - pred).mean()),
                "Spearman": float(spearmanr(y, pred).statistic),
            })
    pl.DataFrame(metrics).sort("split_id", "method", "fold").write_csv(out / "metrics.csv")
    oof.write_parquet(out / "oof.parquet")
    (out / "summary.json").write_text(json.dumps({
        "config": cfg, "rows": oof.height,
        "qualification": "Atomic subset only; no claim about full submitted ensemble",
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
