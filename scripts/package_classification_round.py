"""Package the selected, already audited classification means without refitting."""

import json
import os
from pathlib import Path

import polars as pl

from cyp_challenge.tracks import ENDPOINTS, validate

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
pred = pl.read_parquet(inp / cfg["predictions"])
cuts = pl.read_csv(inp / cfg["thresholds"])
ref = pl.read_csv(inp / cfg["reference"])
aliases = (
    pl.read_parquet(inp / cfg["aliases"])
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
identity = ref.select("Molecule_Name", "SMILES").join(
    aliases, on="Molecule_Name", how="left", validate="1:1"
)
assert identity.height == identity["compound_id"].n_unique() == 750
assert identity["compound_id"].null_count() == 0
parts, summary = [], {}
assert set(cfg["selected"]) == set(ENDPOINTS["classification"])
for ep, method in cfg["selected"].items():
    p = pred.filter((pl.col("endpoint") == ep) & (pl.col("method") == method))
    c = cuts.filter((pl.col("endpoint") == ep) & (pl.col("method") == method))
    assert c["fold"].sort().to_list() == list(range(5))
    threshold = c["threshold"].median()
    assert p.height == p["compound_id"].n_unique() == 750
    assert p["prediction"].is_finite().all() and p["prediction"].is_between(0, 1).all()
    assert (p["threshold"] == threshold).all()
    assert (p["binary"] == (p["prediction"] >= threshold)).all()
    assert set(p["compound_id"]) == set(identity["compound_id"])
    parts.append(p)
    summary[ep] = dict(method=method, threshold=threshold, positives=int(p["binary"].sum()))
long = pl.concat(parts).join(identity, on="compound_id", how="left", validate="m:1")
long.write_parquet(out / "selected-probabilities.parquet")
wide = long.pivot(on="endpoint", index=["SMILES", "Molecule_Name"], values="binary")
wide = wide.select("SMILES", "Molecule_Name", *ENDPOINTS["classification"]).sort("Molecule_Name")
wide.write_parquet(out / "classification.parquet")
report = dict(
    config=cfg,
    endpoints=summary,
    validation=validate(out / "classification.parquet", "classification", inp / cfg["reference"]),
    deployment="Saved fivefold probability means and median inner cutoffs; no additional fit",
)
(out / "summary.json").write_text(json.dumps(report, indent=2))
print(json.dumps(report))
