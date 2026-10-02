"""Audit scored classification probe versus timestamp-matched controls."""

import hashlib
import json
import os
from pathlib import Path

import polars as pl


def row(path):
    frame = pl.read_csv(path).filter(pl.col("Username") == "N283T")
    if frame.height != 1:
        raise ValueError("Expected one participant row")
    return frame.to_dicts()[0]


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    receipt = json.loads((inp / "receipt.json").read_text())
    candidate = pl.read_parquet(inp / "submitted.parquet").sort("Molecule_Name")
    previous = pl.read_parquet(inp / "last.parquet").sort("Molecule_Name")
    best = pl.read_parquet(inp / "best.parquet").sort("Molecule_Name")
    assert receipt["submission_id"] == cfg["submission_id"] and receipt["status"] == "sent"
    assert (
        hashlib.sha256((inp / "submitted.parquet").read_bytes()).hexdigest()
        == receipt["file_sha256"]
    )
    for frame in [candidate, previous, best]:
        assert frame.height == 750 and frame["Molecule_Name"].n_unique() == 750
        assert frame.select("Molecule_Name", "SMILES").equals(
            best.select("Molecule_Name", "SMILES")
        )
    d6, a4 = "CYP2D6_is_TDI", "CYP3A4_is_TDI"
    assert candidate[d6].sum() == 49 and candidate[a4].sum() == 249
    assert candidate[a4].equals(best[a4])
    assert previous[d6].equals(best[d6])
    records, scored = [], {}
    for endpoint in ["MA", "D6", "A4"]:
        after, before = row(inp / f"{endpoint}-after.csv"), row(inp / f"{endpoint}-before.csv")
        assert after["Submitted"] == cfg["expected_scored_timestamp"]
        assert before["Submitted"] == cfg["expected_previous_timestamp"]
        prefix = "MA-" if endpoint == "MA" else ""
        scored[endpoint] = after[prefix + "MCC"]
        for metric in ["MCC", "Accuracy", "Precision", "Recall"]:
            col = prefix + metric
            records.append(
                {
                    "endpoint": endpoint,
                    "metric": metric,
                    "previous": before[col],
                    "current": after[col],
                    "delta": after[col] - before[col],
                    "previous_rank": before["Rank"],
                    "current_rank": after["Rank"],
                }
            )
    best_a4 = row(inp / "A4-best.csv")
    a4_after = row(inp / "A4-after.csv")
    for metric in ["MCC", "Accuracy", "Precision", "Recall"]:
        assert abs(a4_after[metric] - best_a4[metric]) < 1e-12
    assert abs(scored["MA"] - (scored["D6"] + scored["A4"]) / 2) < 1e-12
    metrics = pl.DataFrame(records)
    metrics.write_csv(out / "live-metrics.csv")
    summary = {
        "submission_id": cfg["submission_id"],
        "attempt_id": receipt["id"],
        "sha256": receipt["file_sha256"],
        "scored_timestamp": cfg["expected_scored_timestamp"],
        "scored": True,
        "decision": (
            "Reject D6 50% AP-trunk blend for current Live; retain original D6 performance control"
        ),
        "A4_best_metrics_reproduced": True,
        "macro_arithmetic_checked": True,
        "new_upload": False,
        "metrics": records,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = [
        "# Classification AP-trunk probe: severe D6 Live loss, exact A4 restoration",
        "",
        f"Scored submission `{cfg['submission_id']}` at {cfg['expected_scored_timestamp']}.",
        "",
        "| Endpoint | Metric | Previous | Current | Delta |",
        "|---|---|---:|---:|---:|",
    ]
    for r in records:
        lines.append(
            f"| {r['endpoint']} | {r['metric']} | {r['previous']:.6f} | "
            f"{r['current']:.6f} | {r['delta']:+.6f} |"
        )
    lines.extend(
        [
            "",
            "D6 keeps 49 blind positives but swaps 21 each way. Its ranking AP gain in "
            "conditional labelled-only OOF did not transfer to useful Live calls. The earlier "
            "historical-rate OOF MCC point estimate was slightly negative; it did not resolve "
            "the magnitude of this prospective loss. No hidden molecule labels are inferred.",
            "",
            "A4 is bitwise restored to the prior best 249-call vector and exactly reproduces "
            "its prior MCC/accuracy/precision/recall. This supports the intended restoration "
            "and continuity for A4, not an independent D6 generalization result.",
            "",
            "Reject this D6 composition for current Live. Preserve it as an informative "
            "negative result and retain the original D6 vector as the performance control. "
            "Any rollback is a separate future submission requiring user approval and cooldown. "
            "No new upload, training, extraction or regression change in this audit.",
            "",
        ]
    )
    (out / "report.md").write_text("\n".join(lines))
    print(json.dumps({"scores": scored, "decision": summary["decision"]}, indent=2))


if __name__ == "__main__":
    main()
