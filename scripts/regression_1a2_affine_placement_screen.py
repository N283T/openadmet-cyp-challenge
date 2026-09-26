"""Screen 1A2 affine changes on the current parent bank and a dual-split proxy."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

KEYS = ["compound_id", "endpoint", "fold"]


def score(frame: pl.DataFrame, prediction: np.ndarray) -> dict[str, float]:
    y = frame["y_true"].to_numpy()
    lo = frame["conf_low"].to_numpy()
    hi = frame["conf_high"].to_numpy()
    center = y.mean()
    denominator = (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    return {
        "ST_RAE": float(
            (np.maximum(lo - prediction, 0) + np.maximum(prediction - hi, 0)).sum() / denominator
        ),
        "MAE": float(np.abs(y - prediction).mean()),
        "Spearman": float(spearmanr(y, prediction).statistic),
    }


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    endpoint = cfg["endpoint"]
    full_parent = (
        pl.read_parquet(inp / "parent-oof.parquet")
        .filter((pl.col("endpoint") == endpoint) & (pl.col("method") == "equal_all"))
        .select(KEYS + ["y_true", "conf_low", "conf_high", "prediction"])
        .sort(KEYS)
    )
    proxy = pl.read_parquet(inp / "proxy-oof.parquet").filter(
        (pl.col("endpoint") == endpoint) & (pl.col("method") == cfg["proxy_method"])
    )
    cohorts = [("full_parent", full_parent)]
    for split in ["parent", "scaffold"]:
        part = (
            proxy.filter(pl.col("split_id") == split)
            .select(KEYS + ["y_true", "conf_low", "conf_high", "prediction"])
            .sort(KEYS)
        )
        cohorts.append((f"proxy_{split}", part))
    for name, frame in cohorts:
        if frame.height != 1372 or frame.select(KEYS).n_unique() != frame.height:
            raise ValueError(f"Incomplete {name} OOF")
        if sorted(frame["fold"].unique().to_list()) != list(range(5)):
            raise ValueError(f"Incomplete {name} folds")
    if not full_parent.select("compound_id", "y_true", "conf_low", "conf_high").equals(
        cohorts[1][1].select("compound_id", "y_true", "conf_low", "conf_high")
    ):
        raise ValueError("Full and proxy parent labels differ")
    blind = pl.read_csv(inp / "current.csv")
    if blind.height != 750 or blind["Molecule_Name"].n_unique() != 750:
        raise ValueError("Invalid current submission")
    blind_values = blind[endpoint].to_numpy()

    metrics = []
    fold_metrics = []
    for name, frame in cohorts:
        raw = frame["prediction"].to_numpy()
        center = raw.mean()
        for shift in cfg["shifts"]:
            for scale in cfg["scales"]:
                adjusted = center + shift + scale * (raw - center)
                metrics.append(
                    {
                        "cohort": name,
                        "shift": shift,
                        "scale": scale,
                        "n": frame.height,
                        **score(frame, adjusted),
                    }
                )
                for fold in range(5):
                    mask = frame["fold"].to_numpy() == fold
                    part = frame.filter(pl.Series("keep", mask))
                    fold_metrics.append(
                        {
                            "cohort": name,
                            "fold": fold,
                            "shift": shift,
                            "scale": scale,
                            "n": part.height,
                            **score(part, adjusted[mask]),
                        }
                    )
    pl.DataFrame(metrics).sort("cohort", "shift", "scale").write_csv(out / "metrics.csv")
    pl.DataFrame(fold_metrics).sort("cohort", "fold", "shift", "scale").write_csv(
        out / "fold-metrics.csv"
    )
    pl.concat(
        [frame.with_columns(pl.lit(name).alias("cohort")) for name, frame in cohorts]
    ).write_parquet(out / "base-oof.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "cohorts": {name: frame.height for name, frame in cohorts},
                "blind_mean": float(blind_values.mean()),
                "blind_sample_sd": float(blind_values.std(ddof=1)),
                "qualification": "Only full_parent is the submitted ten-member bank. "
                "The two proxy cohorts average three of its members. "
                "OOF label distribution does not identify blind placement.",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
