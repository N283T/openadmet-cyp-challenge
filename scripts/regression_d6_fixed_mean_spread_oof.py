"""Evaluate the exact proposed D6 spread ratio on keyed held-fold predictions."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def score(frame: pl.DataFrame) -> dict[str, float]:
    y = frame["y_true"].to_numpy()
    p = frame["prediction"].to_numpy()
    lo = frame["conf_low"].to_numpy()
    hi = frame["conf_high"].to_numpy()
    center = y.mean()
    denominator = (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    return {
        "ST_RAE": float((np.maximum(lo - p, 0) + np.maximum(p - hi, 0)).sum() / denominator),
        "MAE": float(np.abs(y - p).mean()),
        "Spearman": float(spearmanr(y, p).statistic),
    }


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    current_path = inp / "submitted.csv"
    variant_path = inp / "variant.csv"
    if sha256(current_path) != cfg["submitted_sha256"]:
        raise ValueError("Current submission SHA mismatch")
    if sha256(variant_path) != cfg["variant_sha256"]:
        raise ValueError("Spread variant SHA mismatch")
    endpoint = cfg["endpoint"]
    current = pl.read_csv(current_path).sort("Molecule_Name")
    variant = pl.read_csv(variant_path).sort("Molecule_Name")
    if current.height != 750 or not current.select("SMILES", "Molecule_Name").equals(
        variant.select("SMILES", "Molecule_Name")
    ):
        raise ValueError("Blind molecule identities differ")
    for column in current.columns[2:]:
        if column != endpoint and not current[column].equals(variant[column]):
            raise ValueError(f"Non-D6 prediction changed: {column}")
    x = current[endpoint].to_numpy()
    z = variant[endpoint].to_numpy()
    ratio = cfg["target_sample_sd"] / x.std(ddof=1)
    expected = (x - x.mean()) * ratio + x.mean()
    if not np.allclose(z, expected, atol=1e-12, rtol=0):
        raise ValueError("Variant is not the expected fixed-mean spread change")

    source = pl.read_parquet(inp / "blend-oof.parquet").filter(
        (pl.col("endpoint") == endpoint)
        & (pl.col("robust_delta") == cfg["robust_delta"])
        & (pl.col("constraint_weight") == cfg["constraint_weight"])
        & (pl.col("blend_weight") == cfg["blend_weight"])
        & (pl.col("placement") == cfg["placement"])
    )
    predictions = []
    folds = []
    for split in ["parent", "scaffold"]:
        for fold in range(5):
            held = source.filter((pl.col("split_id") == split) & (pl.col("fold") == fold))
            if held.height < 200 or held["compound_id"].n_unique() != held.height:
                raise ValueError(f"Incomplete held fold: {split}/{fold}")
            raw = held["prediction"].to_numpy()
            altered = (raw - raw.mean()) * ratio + raw.mean()
            for method, values in [("submitted", raw), ("sd90", altered)]:
                predictions.append(
                    held.select(
                        "compound_id",
                        "endpoint",
                        "split_id",
                        "fold",
                        "y_true",
                        "conf_low",
                        "conf_high",
                    ).with_columns(pl.lit(method).alias("method"), pl.Series("prediction", values))
                )
                folds.append(
                    {
                        "split_id": split,
                        "fold": fold,
                        "method": method,
                        "n": held.height,
                        "mean": float(values.mean()),
                        "sample_sd": float(values.std(ddof=1)),
                        **score(predictions[-1]),
                    }
                )
    oof = pl.concat(predictions).sort("split_id", "method", "compound_id")
    if oof.select("split_id", "method", "compound_id").is_duplicated().any():
        raise ValueError("Duplicate OOF key")
    metrics = []
    for (split, method), part in oof.group_by("split_id", "method"):
        metrics.append({"split_id": split, "method": method, "n": part.height, **score(part)})
    oof.write_parquet(out / "oof.parquet")
    pl.DataFrame(folds).sort("split_id", "fold", "method").write_csv(out / "fold-metrics.csv")
    pl.DataFrame(metrics).sort("split_id", "method").write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "blind_mean": float(x.mean()),
                "blind_original_sample_sd": float(x.std(ddof=1)),
                "blind_variant_sample_sd": float(z.std(ddof=1)),
                "spread_ratio": float(ratio),
                "blind_mean_absolute_change": float(np.abs(z - x).mean()),
                "blind_maximum_absolute_change": float(np.abs(z - x).max()),
                "other_endpoints_exactly_unchanged": True,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
