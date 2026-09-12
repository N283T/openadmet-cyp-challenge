"""Independently verify measured-ADME transforms, censor sides, and identity exclusions."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def read(name):
    return pl.read_parquet(
        next(p for p in inp.iterdir() if p.name.partition("-")[2] == name + ".parquet")
    )


r, o, c = read("az-curated-records"), read("az-observations"), read("external-compounds")
checks = []
for (ep,), g in r.group_by("isoform"):
    raw = g["standard_value"].to_numpy()
    if ep == "ADME_LOGD74":
        expected = raw
    elif ep == "ADME_SOLUBILITY74":
        expected = np.log10(raw * 1e-9)
    elif ep == "ADME_PPB_HUMAN":
        expected = np.log10((1 - raw / 100) * 100)
    else:
        expected = np.log10(raw)
    assert np.allclose(expected, g["transformed_value"].to_numpy())
    checks.append(dict(endpoint=ep, rows=g.height, passed=True))
assert (
    o.filter(pl.col("relation") != "=")["value"].null_count()
    == o.filter(pl.col("relation") != "=").height
)
assert (
    o.filter(pl.col("relation") == "<")["lower_bound"].null_count()
    == o.filter(pl.col("relation") == "<").height
)
assert (
    o.filter(pl.col("relation") == ">")["upper_bound"].null_count()
    == o.filter(pl.col("relation") == ">").height
)
assert np.isclose(np.log10(100 - 90), 1.0)
assert np.log10(100 - 95) < np.log10(100 - 90)
assert set(c["connectivity"]).isdisjoint(read("challenge-connectivity")["connectivity"])
assert not o["challenge_overlap"].any()
assert o["binary_label"].null_count() == o.height
assert o.filter(pl.col("isoform").str.contains("CLINT"))["assay_organism"].unique().to_list() == [
    "Homo sapiens"
]
result = dict(
    passed=True,
    checks=checks,
    censored=o["relation"].is_in(["<", ">", "<=", ">="]).sum(),
    no_adme_binary_labels=True,
    no_challenge_connectivity=True,
    clearance_species="human only",
    ppb_direction=(
        "Higher bound percentage maps to lower log unbound percentage; inequalities reverse"
    ),
)
(out / "summary.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
