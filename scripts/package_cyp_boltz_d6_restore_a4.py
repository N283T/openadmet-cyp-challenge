"""Restore best-scored A4 while preserving the exact reviewed D6 probe."""

import hashlib
import json
import os
from pathlib import Path

import polars as pl

from cyp_challenge.tracks import ENDPOINTS, validate


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    tables = {
        name: pl.read_parquet(inp / f"{name}.parquet").sort("Molecule_Name")
        for name in ["d6-probe", "best-control", "last-sent"]
    }
    reference = (
        pl.read_csv(inp / "reference.csv").select("SMILES", "Molecule_Name").sort("Molecule_Name")
    )
    for name, frame in tables.items():
        validate(inp / f"{name}.parquet", "classification", inp / "reference.csv")
        assert frame.select("SMILES", "Molecule_Name").equals(reference)
    d6, a4 = ENDPOINTS["classification"]
    probe, best, last = [tables[name] for name in ["d6-probe", "best-control", "last-sent"]]
    candidate = probe.with_columns(best[a4])
    assert candidate[d6].equals(probe[d6])
    assert candidate[a4].equals(best[a4])
    assert candidate[d6].sum() == cfg["d6_positive_calls"]
    assert candidate[a4].sum() == cfg["a4_positive_calls"]
    candidate.write_parquet(out / "predictions.parquet")
    candidate.write_csv(out / "predictions.csv")
    validation = validate(out / "predictions.parquet", "classification", inp / "reference.csv")
    changes, changed_frames = [], []
    for name, control in [("last_sent", last), ("best_control", best), ("prior_d6_probe", probe)]:
        for endpoint in ENDPOINTS["classification"]:
            new, old = candidate[endpoint], control[endpoint]
            changes.append(
                {
                    "reference": name,
                    "endpoint": endpoint,
                    "positive_calls": int(new.sum()),
                    "retained": int((new & old).sum()),
                    "added": int((new & ~old).sum()),
                    "removed": int((~new & old).sum()),
                    "changed": int((new != old).sum()),
                }
            )
            changed_frames.append(
                candidate.select("Molecule_Name", "SMILES")
                .with_columns(
                    pl.lit(endpoint).alias("endpoint"),
                    pl.lit(name).alias("reference"),
                    old.alias("previous_call"),
                    new.alias("call"),
                )
                .filter(pl.col("previous_call") != pl.col("call"))
            )
    pl.DataFrame(changes).write_csv(out / "target-changes.csv")
    pl.concat(changed_frames).write_csv(out / "changed-calls.csv")
    summary = {
        "config": cfg,
        "selection_run": cfg["selection_run"],
        "d6_positive_calls": int(candidate[d6].sum()),
        "a4_positive_calls": int(candidate[a4].sum()),
        "D6_bitwise_equal_reviewed_probe": True,
        "A4_bitwise_equal_best_control": True,
        "target_changes": changes,
        "validation": validation,
        "uploaded": False,
        "prediction_sha256": hashlib.sha256((out / "predictions.parquet").read_bytes()).hexdigest(),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = [
        "# Prepared D6 AP probe with exact A4 performance-control restoration",
        "",
        "D6 copies the previously reviewed 50% rank-blend vector exactly: 49 positives, "
        "21 additions and 21 removals versus last sent and best control. No score or rate change.",
        "",
        "A4 copies the best-scored 249-positive vector exactly. Compared with the last sent "
        "275-positive rate probe it removes 26 calls and adds none. No A4 Boltz mixing.",
        "",
        "The fixed 49 D6 calls isolate the score-recipe probe; this is not a claim that "
        "49 is optimal. Both target columns are preserved from immutable source payloads.",
        "",
        "The earlier D6 probe remains retained with its original A4 column. This new package "
        "has a distinct hash and requires its own explicit approval before any upload.",
        "",
        "Local preparation only. No submission, new fitting or extraction, or regression change.",
        "",
    ]
    (out / "report.md").write_text("\n".join(lines))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
