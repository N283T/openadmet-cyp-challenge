"""Audit the reflected Live leaderboard result of the submitted D6-only blend."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def row(path: Path, user: str, submitted: str) -> dict:
    selected = pl.read_csv(path).filter(
        (pl.col("Username") == user) & (pl.col("Submitted") == submitted)
    )
    if selected.height != 1:
        raise ValueError(f"Expected exactly one matched leaderboard row: {path}")
    return selected.to_dicts()[0]


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    receipt = json.loads((inp / "receipt.json").read_text())
    if (
        receipt["submission_id"] != cfg["submission_id"]
        or receipt["id"] != cfg["attempt_id"]
        or receipt["file_sha256"] != cfg["file_sha256"]
        or receipt["status"] != "sent"
    ):
        raise ValueError("Submission receipt mismatch")
    before = pl.read_csv(inp / "before.csv").sort("Molecule_Name")
    after = pl.read_csv(inp / "after.csv").sort("Molecule_Name")
    if before.height != 750 or after.height != 750:
        raise ValueError("Unexpected submission length")
    if not before.select("SMILES", "Molecule_Name").equals(after.select("SMILES", "Molecule_Name")):
        raise ValueError("Molecule identities changed")
    for endpoint in before.columns[2:]:
        if endpoint != cfg["endpoint"] and not before[endpoint].equals(after[endpoint]):
            raise ValueError(f"Unrelated endpoint changed: {endpoint}")
    if sha256(inp / "after.csv") != cfg["file_sha256"]:
        raise ValueError("Submitted payload hash mismatch")
    metrics = []
    saved = {}
    for board in ["d6", "ma"]:
        before_row = row(inp / f"before-{board}.csv", cfg["username"], cfg["before_submitted"])
        after_row = row(inp / f"after-{board}.csv", cfg["username"], cfg["after_submitted"])
        saved[board] = (before_row, after_row)
        for metric in ["ST-RAE", "MAE", "R²", "Spearman's ρ", "Kendall's τ"]:
            column = metric if board == "d6" else "MA-" + metric
            metrics.append(
                {
                    "board": board,
                    "metric": metric,
                    "before": float(before_row[column]),
                    "after": float(after_row[column]),
                    "delta": float(after_row[column] - before_row[column]),
                }
            )
    d6_delta = saved["d6"][1]["ST-RAE"] - saved["d6"][0]["ST-RAE"]
    ma_delta = saved["ma"][1]["MA-ST-RAE"] - saved["ma"][0]["MA-ST-RAE"]
    if not np.isclose(ma_delta, d6_delta / 4, atol=1e-10):
        raise ValueError("MA delta does not match D6-only replacement")
    pl.DataFrame(metrics).write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "rows": before.height,
                "d6_st_rae_delta": d6_delta,
                "ma_st_rae_delta": ma_delta,
                "other_endpoints_exactly_unchanged": True,
                "before_d6_rank": saved["d6"][0]["Rank"],
                "after_d6_rank": saved["d6"][1]["Rank"],
                "before_ma_rank": saved["ma"][0]["Rank"],
                "after_ma_rank": saved["ma"][1]["Rank"],
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
