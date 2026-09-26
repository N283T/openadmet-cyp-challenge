"""Package a 1A2-only narrowing and its combination with prepared D6 SD-0.90."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.tracks import validate


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    current_path = inp / "current.csv"
    d6_path = inp / "d6-sd90.csv"
    if sha256(current_path) != cfg["source_sha256"]:
        raise ValueError("Current scored submission SHA mismatch")
    if sha256(d6_path) != cfg["d6_sd90_sha256"]:
        raise ValueError("D6 spread submission SHA mismatch")
    current = pl.read_csv(current_path).sort("Molecule_Name")
    d6 = pl.read_csv(d6_path).sort("Molecule_Name")
    if current.height != 750 or not current.select("SMILES", "Molecule_Name").equals(
        d6.select("SMILES", "Molecule_Name")
    ):
        raise ValueError("Submission identities mismatch")
    endpoint = cfg["endpoint"]
    for column in current.columns[2:]:
        if column != "CYP2D6_pIC50_direct_inhibition" and not current[column].equals(d6[column]):
            raise ValueError(f"Prepared D6 payload changed another endpoint: {column}")
    original = current[endpoint].to_numpy()
    adjusted = original.mean() + cfg["scale"] * (original - original.mean())
    if adjusted.min() < 1 or not np.isfinite(adjusted).all():
        raise ValueError("Invalid 1A2 adjustment")
    a1_only = current.with_columns(pl.Series(endpoint, adjusted))
    combined = d6.with_columns(pl.Series(endpoint, adjusted))
    for name, frame in [("1a2-only.csv", a1_only), ("1a2-d6-combined.csv", combined)]:
        frame.write_csv(out / name)
        report = validate(out / name, "regression", inp / "reference.csv")
        if not report["valid"]:
            raise ValueError(report)
    a1_only.select("Molecule_Name", "SMILES", endpoint).with_columns(
        pl.lit("blind").alias("split_id"),
        pl.lit(-1).alias("fold"),
        pl.Series("scored_prediction", original),
    ).write_parquet(out / "1a2-comparison.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "rows": current.height,
                "original_mean": float(original.mean()),
                "original_sample_sd": float(original.std(ddof=1)),
                "adjusted_mean": float(adjusted.mean()),
                "adjusted_sample_sd": float(adjusted.std(ddof=1)),
                "mean_absolute_change": float(np.abs(adjusted - original).mean()),
                "maximum_absolute_change": float(np.abs(adjusted - original).max()),
                "a1_only_sha256": sha256(out / "1a2-only.csv"),
                "combined_sha256": sha256(out / "1a2-d6-combined.csv"),
                "other_endpoints_unchanged_in_a1_only": True,
                "combined_changes_only_a1_and_d6": True,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
