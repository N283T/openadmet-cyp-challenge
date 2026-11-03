"""Re-threshold the deployed L4 rank-mean D6 score at another call count; A4 stays unchanged.

The score is rebuilt exactly as in `package_cyp_d6_gate_candidate.py` (equal percentile-rank mean of
the deployed Sep 15 blind score and the frozen `pool_lgbm` blind score). The top-49 calls must
reproduce the D6 column of the last sent payload, so only the number of calls changes. Nothing is
fitted or uploaded.
"""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from package_cyp_d6_gate_candidate import (
    IDENTIFIERS,
    blind_column,
    payload,
    percentile_rank,
    top_calls,
)

from cyp_challenge.tracks import validate


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    d6, a4 = cfg["d6_endpoint"], cfg["a4_endpoint"]
    reference = pl.read_csv(inp / "reference.csv").select(IDENTIFIERS).sort("Molecule_Name")
    names = (
        pl.read_parquet(inp / "test_mapping.parquet")
        .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
        .unique()
    )
    last = payload(inp / "last_sent.parquet", reference)
    assert int(last[d6].sum()) == cfg["reference_calls"]
    assert int(last[a4].sum()) == cfg["a4_positive_calls"]

    legacy = pl.read_parquet(inp / "sep15_blind_scores.parquet").filter(
        pl.col("model") == cfg["legacy_model"]
    )
    assert legacy.height == 750 and legacy["score"].is_finite().all()
    table = legacy.select("compound_id", pl.col("score").alias("legacy_score")).join(
        names, on="compound_id", how="left", validate="1:1"
    )
    table = table.join(
        blind_column(pl.read_parquet(inp / "test_scores.parquet"), cfg["pool_column"]),
        on="compound_id",
        validate="1:1",
    ).sort("compound_id")
    score = 0.5 * (
        percentile_rank(table["legacy_score"].to_numpy())
        + percentile_rank(table[cfg["pool_column"]].to_numpy())
    )
    last_call = (
        table.select("Molecule_Name")
        .join(last.select("Molecule_Name", pl.col(d6)), on="Molecule_Name", maintain_order="left")[
            d6
        ]
        .to_numpy()
    )
    assert np.array_equal(top_calls(score, cfg["reference_calls"]), last_call), (
        "score not reproduced"
    )

    count = cfg["d6_positive_calls"]
    calls = top_calls(score, count)
    order = np.argsort(-score, kind="stable")
    assert score[order[count - 1]] - score[order[count]] > 0, "Tied call boundary"
    assert (calls | ~last_call).all(), "new calls must contain the reference calls"

    d6_calls = table.select("Molecule_Name").with_columns(pl.Series(d6, calls))
    candidate = (
        last.drop(d6)
        .join(d6_calls, on="Molecule_Name", validate="1:1", maintain_order="left")
        .select(*IDENTIFIERS, d6, a4)
    )
    assert candidate[a4].equals(last[a4]) and int(candidate[d6].sum()) == count
    candidate.write_parquet(out / "predictions.parquet")
    candidate.write_csv(out / "predictions.csv")
    validation = validate(out / "predictions.parquet", "classification", inp / "reference.csv")
    summary = {
        "d6_calls": count,
        "reference_calls": cfg["reference_calls"],
        "added_calls": int(calls.sum() - last_call.sum()),
        "a4_unchanged": True,
        "boundary_gap": float(score[order[count - 1]] - score[order[count]]),
        "validation": validation,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps(summary, default=str))


if __name__ == "__main__":
    main()
