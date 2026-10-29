"""Take selected columns from a new payload and keep every other column from a base payload."""

import json
import os
from pathlib import Path

import polars as pl


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    base = pl.read_csv(inp / "base.csv")
    new = pl.read_csv(inp / "new.csv")
    assert base["Molecule_Name"].equals(new["Molecule_Name"]) and base["SMILES"].equals(new["SMILES"])
    payload = base
    for column in cfg["take_from_new"]:
        payload = payload.with_columns(new[column])
    payload.write_csv(out / "predictions.csv")
    summary = {
        column: {
            "changed_rows": int((base[column] != payload[column]).sum()),
            "mean_abs": float((base[column] - payload[column]).abs().mean()),
        }
        for column in base.columns
        if column.endswith("_direct_inhibition")
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
