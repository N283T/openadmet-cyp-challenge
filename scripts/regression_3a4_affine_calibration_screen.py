"""Screen fixed location and spread changes to the scored CYP3A4 OOF."""

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
    frame = (
        pl.read_parquet(inp / "scored-bank-oof.parquet")
        .filter((pl.col("endpoint") == cfg["endpoint"]) & (pl.col("method") == cfg["oof_method"]))
        .sort(["compound_id", "endpoint", "fold"])
    )
    if frame.height != 2300 or frame.select("compound_id", "endpoint", "fold").n_unique() != 2300:
        raise ValueError("Incomplete scored 3A4 OOF")
    blind = pl.read_csv(inp / "current.csv")
    if blind.height != 750 or blind["Molecule_Name"].n_unique() != 750:
        raise ValueError("Invalid current submission")
    raw = frame["prediction"].to_numpy()
    blind_raw = blind[cfg["endpoint"]].to_numpy()
    folds = frame["fold"].to_numpy()
    rows = []
    fold_rows = []
    for shift in cfg["shifts"]:
        for scale in cfg["scales"]:
            adjusted = raw.mean() + shift + scale * (raw - raw.mean())
            rows.append(
                {"shift": shift, "scale": scale, "n": frame.height, **score(frame, adjusted)}
            )
            for fold in sorted(np.unique(folds)):
                mask = folds == fold
                part = frame.filter(pl.Series("keep", mask))
                fold_rows.append(
                    {
                        "shift": shift,
                        "scale": scale,
                        "fold": int(fold),
                        "n": part.height,
                        **score(part, adjusted[mask]),
                    }
                )
    pl.DataFrame(rows).sort("ST_RAE").write_csv(out / "metrics.csv")
    pl.DataFrame(fold_rows).sort("shift", "scale", "fold").write_csv(out / "fold-metrics.csv")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "oof_rows": frame.height,
                "oof_mean": float(raw.mean()),
                "oof_sample_sd": float(raw.std(ddof=1)),
                "blind_mean": float(blind_raw.mean()),
                "blind_sample_sd": float(blind_raw.std(ddof=1)),
                "qualification": "OOF screen is developmental; blind placement may differ.",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
