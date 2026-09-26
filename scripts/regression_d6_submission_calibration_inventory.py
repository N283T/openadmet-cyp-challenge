"""Align scored D6 submission vectors and identify pure affine comparisons."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    endpoint = cfg["endpoint"]
    frames = {}
    vectors = {}
    rows = []
    keyed = []
    for item in cfg["submissions"]:
        alias = item["alias"]
        path = inp / alias
        frame = (pl.read_parquet(path) if path.suffix == ".parquet" else pl.read_csv(path)).sort(
            "Molecule_Name"
        )
        if frame.height != 750 or frame["Molecule_Name"].n_unique() != 750:
            raise ValueError(f"Invalid molecule coverage: {alias}")
        if frames and not frame.select("Molecule_Name", "SMILES").equals(
            next(iter(frames.values())).select("Molecule_Name", "SMILES")
        ):
            raise ValueError(f"Molecule mismatch: {alias}")
        frames[alias] = frame
        y = frame[endpoint].to_numpy()
        if not np.isfinite(y).all():
            raise ValueError(f"Nonfinite prediction: {alias}")
        vectors[alias] = y
        rows.append(
            {
                "alias": alias,
                "submission_id": item["id"],
                "mean": float(y.mean()),
                "sd_sample": float(y.std(ddof=1)),
                "minimum": float(y.min()),
                "maximum": float(y.max()),
            }
        )
        keyed.append(
            frame.select("Molecule_Name", "SMILES").with_columns(
                pl.lit(item["id"]).alias("submission_id"),
                pl.lit(endpoint).alias("endpoint"),
                pl.lit("blind").alias("split_id"),
                pl.lit(-1).alias("fold"),
                pl.Series("prediction", y.astype(np.float64)),
            )
        )
    comparisons = []
    aliases = list(vectors)
    for i, left in enumerate(aliases):
        x = vectors[left]
        for right in aliases[i + 1 :]:
            y = vectors[right]
            scale = float(np.cov(x, y, ddof=0)[0, 1] / np.var(x))
            offset = float(y.mean() - scale * x.mean())
            residual = y - (scale * x + offset)
            comparisons.append(
                {
                    "left": left,
                    "right": right,
                    "pearson": float(np.corrcoef(x, y)[0, 1]),
                    "spearman": float(spearmanr(x, y).statistic),
                    "affine_scale": scale,
                    "affine_offset": offset,
                    "affine_rmse": float(np.sqrt(np.mean(residual**2))),
                    "affine_max_error": float(np.abs(residual).max()),
                    "exact_affine_1e8": bool(np.abs(residual).max() < 1e-8),
                }
            )
    pl.DataFrame(rows).write_csv(out / "moments.csv")
    pl.DataFrame(comparisons).write_csv(out / "affine-comparisons.csv")
    pl.concat(keyed).write_parquet(out / "aligned-predictions.parquet")
    (out / "summary.json").write_text(json.dumps({"config": cfg}, indent=2) + "\n")


if __name__ == "__main__":
    main()
