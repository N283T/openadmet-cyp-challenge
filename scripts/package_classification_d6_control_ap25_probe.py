"""Package a four-swap D6 control/AP probe with the validated A4 calls."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.tracks import validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
base = pl.read_parquet(inp / cfg["base_submission"]).sort("Molecule_Name")
scores = pl.read_parquet(inp / cfg["decision_test"])
aliases = (
    pl.read_parquet(inp / cfg["aliases"])
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
parts = [
    scores.filter(pl.col("model") == name).sort("compound_id")
    for name in [cfg["control_model"], cfg["ap_model"]]
]
if any(part.height != 750 for part in parts) or not parts[0]["compound_id"].equals(
    parts[1]["compound_id"]
):
    raise ValueError("Blind model scores are incomplete or misaligned")
if (
    not np.isfinite(parts[0]["score"].to_numpy()).all()
    or not np.isfinite(parts[1]["score"].to_numpy()).all()
):
    raise ValueError("Nonfinite blind score")
if int(parts[0]["call"].sum()) != cfg["positive_calls"]:
    raise ValueError("Control positive-call count changed")
blended = (1 - cfg["ap_weight"]) * parts[0]["score"].to_numpy() + cfg["ap_weight"] * parts[1][
    "score"
].to_numpy()
chosen = np.zeros(750, dtype=bool)
chosen[np.argsort(-blended, kind="stable")[: cfg["positive_calls"]]] = True
decision = (
    parts[0]
    .select("compound_id")
    .with_columns(
        pl.lit("sep15_ap25_rank_blend").alias("model"),
        pl.Series("score", blended),
        pl.Series("call", chosen),
    )
)
decision.write_parquet(out / "decision-test.parquet")
mapped = aliases.join(
    decision.select("compound_id", pl.col("call").alias("new_d6")),
    on="compound_id",
    validate="1:1",
)
if mapped.height != 750 or mapped["Molecule_Name"].n_unique() != 750:
    raise ValueError("Blind package mapping incomplete")
candidate = (
    base.join(mapped.select("Molecule_Name", "new_d6"), on="Molecule_Name", validate="1:1")
    .with_columns(pl.col("new_d6").alias("CYP2D6_is_TDI"))
    .drop("new_d6")
    .select(base.columns)
    .sort("Molecule_Name")
)
if not candidate["CYP3A4_is_TDI"].equals(base["CYP3A4_is_TDI"]):
    raise ValueError("CYP3A4 changed")
old = base["CYP2D6_is_TDI"].to_numpy().astype(bool)
new = candidate["CYP2D6_is_TDI"].to_numpy().astype(bool)
swaps = {
    "removed": int((old & ~new).sum()),
    "added": int((~old & new).sum()),
    "retained": int((old & new).sum()),
}
candidate.write_parquet(out / "predictions.parquet")
candidate.write_csv(out / "predictions.csv")
validation = validate(out / "predictions.parquet", "classification", inp / cfg["reference"])
summary = {
    "model": "sep15_ap25_rank_blend",
    "d6_positive_calls": int(new.sum()),
    "d6_swaps": swaps,
    "a4_positive_calls": int(candidate["CYP3A4_is_TDI"].sum()),
    "a4_bitwise_unchanged": True,
    "validation": validation,
    "prediction_sha256": hashlib.sha256((out / "predictions.parquet").read_bytes()).hexdigest(),
    "uploaded": False,
    "evidence_gate": "INSUFFICIENT_EVIDENCE",
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
