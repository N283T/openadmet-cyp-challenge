"""Package fixed placement controls with official identity/schema validation."""

import hashlib
import json
import os
from pathlib import Path

import polars as pl

from cyp_challenge.tracks import ENDPOINTS, validate

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


scenarios = pl.read_parquet(path(0))
aliases = (
    pl.read_parquet(path(1))
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
reference = pl.read_csv(path(2))
identity = reference.select("Molecule_Name", "SMILES").join(
    aliases, on="Molecule_Name", validate="1:1"
)
assert len(identity) == identity["compound_id"].n_unique() == 750
expected = (
    identity.select("compound_id")
    .join(pl.DataFrame({"endpoint": ENDPOINTS["regression"]}), how="cross")
    .sort("compound_id", "endpoint")
)
manifest = []
for name in ["raw", "author_solved_strae", "d6_only_author_strae"]:
    h = (
        scenarios.filter(pl.col("scenario") == name)
        .select("compound_id", "endpoint", "prediction")
        .sort("compound_id", "endpoint")
    )
    assert h.select("compound_id", "endpoint").equals(expected)
    assert h["prediction"].is_finite().all()
    wide = (
        h.join(identity, on="compound_id", validate="m:1")
        .pivot(on="endpoint", index=["Molecule_Name", "SMILES"], values="prediction")
        .select("SMILES", "Molecule_Name", *ENDPOINTS["regression"])
        .sort("Molecule_Name")
    )
    filename = out / f"{name}.csv"
    wide.write_csv(filename)
    checked = validate(filename, "regression", path(2))
    assert checked["valid"], checked
    manifest.append(
        dict(
            scenario=name,
            filename=filename.name,
            rows=len(wide),
            sha256=hashlib.sha256(filename.read_bytes()).hexdigest(),
            validation=checked,
        )
    )
(out / "manifest.json").write_text(
    json.dumps(
        dict(
            files=manifest,
            source="same previously submitted vector; source assumption calibration; no upload",
        ),
        indent=2,
    )
)
