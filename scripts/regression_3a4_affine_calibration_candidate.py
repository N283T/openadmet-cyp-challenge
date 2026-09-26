"""Apply a fixed 3A4 affine map without changing molecular ordering."""

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
    if sha256(inp / "current.csv") != cfg["source_sha256"]:
        raise ValueError("Scored control SHA mismatch")
    if sha256(inp / "combined.csv") != cfg["combined_sha256"]:
        raise ValueError("Prepared combined SHA mismatch")
    current = pl.read_csv(inp / "current.csv").sort("Molecule_Name")
    combined = pl.read_csv(inp / "combined.csv").sort("Molecule_Name")
    if current.height != 750 or not current.select("Molecule_Name", "SMILES").equals(
        combined.select("Molecule_Name", "SMILES")
    ):
        raise ValueError("Submission identities mismatch")
    endpoint = cfg["endpoint"]
    if not current[endpoint].equals(combined[endpoint]):
        raise ValueError("Prepared combination already changed 3A4")
    original = current[endpoint].to_numpy()
    adjusted = original.mean() + cfg["shift"] + cfg["scale"] * (original - original.mean())
    if adjusted.min() < 1 or not np.isfinite(adjusted).all():
        raise ValueError("Invalid corrected 3A4 range")
    outputs = {
        "3a4-only.csv": current.with_columns(pl.Series(endpoint, adjusted)),
        "1a2-d6-3a4-combined.csv": combined.with_columns(pl.Series(endpoint, adjusted)),
    }
    for name, frame in outputs.items():
        frame.write_csv(out / name)
        report = validate(out / name, "regression", inp / "reference.csv")
        if not report["valid"]:
            raise ValueError(report)
    only = outputs["3a4-only.csv"]
    for column in current.columns[2:]:
        if column != endpoint and not only[column].equals(current[column]):
            raise ValueError(f"3A4-only changed another endpoint: {column}")
    for column in combined.columns[2:]:
        if column != endpoint and not outputs["1a2-d6-3a4-combined.csv"][column].equals(
            combined[column]
        ):
            raise ValueError(f"Triple combination changed another endpoint: {column}")
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
                "only_sha256": sha256(out / "3a4-only.csv"),
                "combined_sha256": sha256(out / "1a2-d6-3a4-combined.csv"),
                "other_endpoints_unchanged_in_only": True,
                "existing_1a2_d6_changes_preserved_in_combined": True,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
