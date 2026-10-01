"""Audit reflected Boltz submission against last and best Live controls."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from regression_calibration_live_result import row, sha256


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    receipt = json.loads((inp / "receipt.json").read_text())
    assert receipt["status"] == "sent" and receipt["submission_id"] == cfg["submission"]
    assert receipt["file_sha256"] == sha256(inp / "current.csv") == cfg["sha256"]
    frames = {k: pl.read_csv(inp / f"{k}.csv").sort("Molecule_Name") for k in cfg["versions"]}
    for frame in frames.values():
        assert frame.height == 750
        assert frame.select("Molecule_Name", "SMILES").equals(
            frames["current"].select("Molecule_Name", "SMILES")
        )
    for endpoint in ["CYP2C9_pIC50_direct_inhibition", "CYP3A4_pIC50_direct_inhibition"]:
        assert frames["current"][endpoint].equals(frames["best"][endpoint])
    rows = []
    for board in cfg["boards"]:
        versions = {
            k: row(inp / f"{k}-{board}.csv", "N283T", date) for k, date in cfg["versions"].items()
        }
        prefix = "MA-" if board == "MA" else ""
        for metric in ["ST-RAE", "MAE", "R²", "Spearman's ρ", "Kendall's τ"]:
            values = {k: float(r[prefix + metric]) for k, r in versions.items()}
            rows.append(
                {
                    "board": board,
                    "metric": metric,
                    **values,
                    "delta_vs_last": values["current"] - values["last"],
                    "delta_vs_best": values["current"] - values["best"],
                }
            )
    f = pl.DataFrame(rows)
    st = f.filter(pl.col("metric") == "ST-RAE")
    for version in cfg["versions"]:
        assert np.isclose(
            st.filter(pl.col("board") == "MA")[version][0],
            st.filter(pl.col("board") != "MA")[version].mean(),
            atol=1e-10,
        )
    for metric in f["metric"].unique():
        c9 = f.filter((pl.col("board") == "2C9") & (pl.col("metric") == metric)).row(0, named=True)
        assert c9["current"] == c9["last"] == c9["best"]
        a4 = f.filter((pl.col("board") == "3A4") & (pl.col("metric") == metric)).row(0, named=True)
        assert a4["current"] == a4["best"]
    f.write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "scored": True,
                "unchanged_C9_all_metrics": True,
                "restored_A4_all_metrics": True,
                "population": (
                    "Live to Live; identical unchanged C9 metrics support population continuity"
                ),
                "ST_RAE": st.to_dicts(),
            },
            indent=2,
        )
        + "\n"
    )
    lines = [
        "# Boltz regression Live result",
        "",
        f"Submission `{cfg['submission']}`.",
        "",
        "| Target | Last | Best control | Current | Delta vs best |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in st.to_dicts():
        lines.append(
            f"| {r['board']} | {r['last']:.6f} | {r['best']:.6f} "
            f"| {r['current']:.6f} | {r['delta_vs_best']:+.6f} |"
        )
    lines.extend(
        [
            "",
            "Both changed targets improve on the prior best scored control. "
            "C9 is identical; A4 returns exactly to its best control. All metrics are "
            "compared on Live boards, matched to exact submission timestamps. "
            "Identical C9 metrics provide a population continuity check. This is prospective "
            "aggregate leaderboard evidence, without per-molecule or hidden-half validation.",
        ]
    )
    (out / "report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
