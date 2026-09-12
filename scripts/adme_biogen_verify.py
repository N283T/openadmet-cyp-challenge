"""Verify Biogen logged values are unchanged, provenance resolves, and exclusions hold."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def path(name):
    return next(p for p in inp.iterdir() if p.name.partition("-")[2] == name)


def read(name):
    return pl.read_parquet(path(name + ".parquet"))


o = read("biogen-observations")
c = read("external-compounds")
assert np.array_equal(o["value"].to_numpy(), o["raw_value"].cast(pl.Float64).to_numpy())
assert o["relation"].unique().to_list() == ["="]
assert o["binary_label"].null_count() == o.height
assert not o["challenge_overlap"].any()
assert set(c["connectivity"]).isdisjoint(read("challenge-connectivity")["connectivity"])
assert o["source_csv_sha256"].unique().to_list() == [
    hashlib.sha256(path("ADME_public_set_3521.csv").read_bytes()).hexdigest()
]
assert o["assay_id"].n_unique() == 6
assert o.select("compound_id", "assay_id").n_unique() == o.height
result = dict(
    passed=True,
    unchanged_logged_values=o.height,
    heads=6,
    no_double_log=True,
    no_fabricated_binary=True,
    no_challenge_connectivity=True,
    source_hash_verified=True,
    no_explicit_censors=True,
)
(out / "summary.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
