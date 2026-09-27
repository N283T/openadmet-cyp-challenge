"""Compare the accepted regression calibration probe to its scored control."""

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
        raise ValueError(f"Expected one matched leaderboard row: {path}")
    return selected.to_dicts()[0]


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    receipt = json.loads((inp / "receipt.json").read_text())
    if any(
        receipt[key] != expected
        for key, expected in [
            ("submission_id", cfg["submission_id"]),
            ("id", cfg["attempt_id"]),
            ("file_sha256", cfg["file_sha256"]),
            ("status", "sent"),
        ]
    ):
        raise ValueError("Submission receipt mismatch")
    if sha256(inp / "after.csv") != cfg["file_sha256"]:
        raise ValueError("Submitted payload hash mismatch")
    before = pl.read_csv(inp / "before.csv").sort("Molecule_Name")
    after = pl.read_csv(inp / "after.csv").sort("Molecule_Name")
    if (
        before.height != 750
        or after.height != 750
        or not before.select("Molecule_Name", "SMILES").equals(
            after.select("Molecule_Name", "SMILES")
        )
    ):
        raise ValueError("Molecule identities differ")
    changes = {
        endpoint: int(np.count_nonzero(before[endpoint].to_numpy() != after[endpoint].to_numpy()))
        for endpoint in before.columns[2:]
    }
    if changes != {
        "CYP1A2_pIC50_direct_inhibition": 750,
        "CYP2C9_pIC50_direct_inhibition": 0,
        "CYP2D6_pIC50_direct_inhibition": 750,
        "CYP3A4_pIC50_direct_inhibition": 750,
    }:
        raise ValueError(f"Unexpected endpoint changes: {changes}")
    metrics = []
    summary = {"config": cfg, "target_changes": changes, "boards": {}}
    deltas = {}
    for board in cfg["boards"]:
        old = row(inp / f"before-{board}.csv", cfg["username"], cfg["before_submitted"])
        new = row(inp / f"after-{board}.csv", cfg["username"], cfg["after_submitted"])
        prefix = "MA-" if board == "MA" else ""
        deltas[board] = new[prefix + "ST-RAE"] - old[prefix + "ST-RAE"]
        summary["boards"][board] = {
            "before_rank": old["Rank"],
            "after_rank": new["Rank"],
            "before_st_rae": old[prefix + "ST-RAE"],
            "after_st_rae": new[prefix + "ST-RAE"],
            "delta_st_rae": deltas[board],
        }
        for metric in ["ST-RAE", "MAE", "R²", "Spearman's ρ", "Kendall's τ"]:
            column = prefix + metric
            metrics.append(
                {
                    "board": board,
                    "metric": metric,
                    "before": float(old[column]),
                    "after": float(new[column]),
                    "delta": float(new[column] - old[column]),
                }
            )
    if not np.isclose(deltas["MA"], sum(deltas[x] for x in cfg["boards"][1:]) / 4, atol=1e-10):
        raise ValueError("Macro delta does not equal mean endpoint delta")
    if deltas["2C9"] != 0:
        raise ValueError("Unchanged 2C9 endpoint score moved")
    pl.DataFrame(metrics).write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
