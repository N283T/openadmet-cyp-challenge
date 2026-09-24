"""Package AP-mass90 D6 calls while preserving the scored A4 vector."""

import hashlib
import json
import os
from pathlib import Path

import polars as pl

from cyp_challenge.tracks import validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
D6, A4 = "CYP2D6_is_TDI", "CYP3A4_is_TDI"


def read(key: str) -> pl.DataFrame:
    path = inp / cfg[key]
    return pl.read_csv(path) if path.suffix == ".csv" else pl.read_parquet(path)


base = read("base_submission").sort("Molecule_Name")
aliases = (
    read("aliases")
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
selected = (
    read("decision_test")
    .filter(pl.col("model") == cfg["d6_model"])
    .select("compound_id", pl.col("call").cast(pl.Boolean).alias("new_d6"))
)
replacement = aliases.join(selected, on="compound_id", validate="1:1").select(
    "Molecule_Name", "new_d6"
)
candidate = (
    base.join(replacement, on="Molecule_Name", validate="1:1")
    .with_columns(pl.col("new_d6").alias(D6))
    .drop("new_d6")
    .select(base.columns)
    .sort("Molecule_Name")
)
if not candidate[A4].equals(base[A4]):
    raise ValueError("CYP3A4 changed")
if int(candidate[D6].sum()) != cfg["d6_positive_calls"]:
    raise ValueError("Unexpected CYP2D6 call count")
old_d6 = base[D6].to_numpy().astype(bool)
new_d6 = candidate[D6].to_numpy().astype(bool)
swaps = {
    "base_only": int((old_d6 & ~new_d6).sum()),
    "candidate_only": int((~old_d6 & new_d6).sum()),
    "both_positive": int((old_d6 & new_d6).sum()),
}
candidate.write_parquet(out / "predictions.parquet")
candidate.write_csv(out / "predictions.csv")
validation = validate(out / "predictions.parquet", "classification", inp / cfg["reference"])
sha256 = hashlib.sha256((out / "predictions.parquet").read_bytes()).hexdigest()
summary = {
    "D6_model": cfg["d6_model"],
    "D6_positive_calls": int(candidate[D6].sum()),
    "D6_swaps": swaps,
    "A4_positive_calls": int(candidate[A4].sum()),
    "A4_bitwise_unchanged": True,
    "validation": validation,
    "prediction_sha256": sha256,
    "uploaded": False,
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
