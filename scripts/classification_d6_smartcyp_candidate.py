"""Package the audited D6 SMARTCyp blend while restoring current A4 exactly."""

import json
import os
from pathlib import Path

import polars as pl

from cyp_challenge.tracks import ENDPOINTS, validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
d6 = "CYP2D6_is_TDI"
a4 = "CYP3A4_is_TDI"

mechanism = pl.read_parquet(inp / cfg["mechanism_candidate"]).sort("Molecule_Name")
current = pl.read_parquet(inp / cfg["current_candidate"]).sort("Molecule_Name")
if not mechanism.select("SMILES", "Molecule_Name").equals(
    current.select("SMILES", "Molecule_Name")
):
    raise ValueError("Candidate identities differ")

candidate = mechanism.drop(a4).with_columns(current[a4]).select(
    "SMILES", "Molecule_Name", *ENDPOINTS["classification"]
)
if not candidate[a4].equals(current[a4]):
    raise ValueError("A4 was not restored bitwise")

old_d6 = current[d6].to_numpy().astype(bool)
new_d6 = candidate[d6].to_numpy().astype(bool)
candidate.write_parquet(out / "predictions.parquet")
candidate.write_csv(out / "predictions.csv")
summary = {
    "source_run": cfg["source_run"],
    "selection": cfg["selection"],
    "positive_calls": {d6: int(candidate[d6].sum()), a4: int(candidate[a4].sum())},
    "d6_current_only": int((old_d6 & ~new_d6).sum()),
    "d6_candidate_only": int((~old_d6 & new_d6).sum()),
    "a4_bitwise_unchanged": True,
    "validation": validate(
        out / "predictions.parquet", "classification", inp / cfg["reference"]
    ),
    "uploaded": False,
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
