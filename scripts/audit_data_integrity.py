"""Independent raw-data freshness, molecular identity, label and submission audit."""

import hashlib
import json
import os
from pathlib import Path

import httpx
import numpy as np
import polars as pl
from rdkit import Chem

from cyp_challenge.preprocessing import standardize

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])


def path(name):
    return next(inp.glob("*-" + name))


def read(name):
    return pl.read_parquet(path(name))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


summary = {}
with httpx.Client(timeout=120, follow_redirects=True) as client:
    response = client.get("https://huggingface.co/api/datasets/openadmet/cyp-challenge-train-test")
    response.raise_for_status()
    meta = response.json()
    (out / "hf-metadata.json").write_text(json.dumps(meta, indent=2))
    revision = meta["sha"]
    fresh = []
    for local in sorted(inp.glob("*-*.csv")):
        name = local.name.split("-", 1)[1]
        response = client.get(
            f"https://huggingface.co/datasets/openadmet/cyp-challenge-train-test/resolve/{revision}/{name}"
        )
        response.raise_for_status()
        remote = out / name
        remote.write_bytes(response.content)
        a, b = (
            pl.read_csv(local, infer_schema_length=None),
            pl.read_csv(remote, infer_schema_length=None),
        )
        fresh.append(
            dict(
                file=name,
                local_sha=sha(local),
                remote_sha=sha(remote),
                exact_bytes=sha(local) == sha(remote),
                local_rows=a.height,
                remote_rows=b.height,
                table_equal=a.equals(b),
            )
        )
    summary["freshness"] = dict(revision=revision, files=fresh)
    print("FRESHNESS", json.dumps(summary["freshness"]), flush=True)

aliases, compounds, saved = (
    read("aliases.parquet"),
    read("compounds.parquet"),
    read("standardization.parquet"),
)
recomputed = []
for i, s in enumerate(saved["raw_smiles"]):
    r = standardize(s)
    raw = Chem.MolFromSmiles(s)
    parent = Chem.MolFromSmiles(r["standardized_smiles"]) if r["standardized_smiles"] else None
    r["raw_heavy_atoms"] = raw.GetNumHeavyAtoms() if raw else None
    r["parent_heavy_atoms"] = parent.GetNumHeavyAtoms() if parent else None
    r["raw_fragments"] = len(Chem.GetMolFrags(raw)) if raw else None
    recomputed.append(r)
    if i % 1000 == 0:
        print("SMILES", i, flush=True)
r = pl.DataFrame(recomputed)
r.write_parquet(out / "recomputed-standardization.parquet")
for col in ["compound_id", "standardized_smiles", "status"]:
    assert saved.sort("raw_smiles")[col].equals(r.sort("raw_smiles")[col]), col
valid = r.filter(pl.col("compound_id").is_not_null())
assert compounds.sort("compound_id").equals(
    valid.select("compound_id", "standardized_smiles").unique().sort("compound_id")
)
changed = r.filter(pl.col("status") == "changed")
changed.write_csv(out / "changed-smiles.csv")
name_conflicts = (
    aliases.group_by("molecule_name")
    .agg(pl.col("compound_id").n_unique().alias("n"))
    .filter(pl.col("n") > 1)
)
name_conflicts.write_csv(out / "name-conflicts.csv")
summary["smiles"] = dict(
    unique_raw=r.height,
    compounds=compounds.height,
    statuses=r.group_by("status").len().to_dicts(),
    name_conflicts=name_conflicts.height,
    recomputation_exact=True,
    heavy_atom_changes=changed.filter(
        pl.col("raw_heavy_atoms") != pl.col("parent_heavy_atoms")
    ).height,
)

measurements = read("measurements.parquet")
raw_labels = []
for name in [
    "cyp-challenge-TRAIN_inhibition.csv",
    "cyp-challenge-TRAIN_TDI.csv",
    "cyp-challenge-TRAIN_Emax.csv",
]:
    frame = pl.read_csv(path(name), infer_schema_length=None)
    mapping = (
        aliases.filter(pl.col("source_file") == name)
        .select(
            pl.col("molecule_name").alias("Molecule_Name"),
            pl.col("raw_smiles").alias("SMILES"),
            "compound_id",
        )
        .unique()
    )
    frame = frame.join(mapping, on=["Molecule_Name", "SMILES"], how="left", validate="m:1")
    for ep in [c for c in frame.columns if c.endswith("_pIC50_direct_inhibition")]:
        raw_labels.append(
            frame.filter(pl.col(ep).is_not_null() & pl.col("compound_id").is_not_null()).select(
                "compound_id",
                pl.lit(ep).alias("endpoint"),
                pl.col(ep).alias("value"),
                pl.col(ep + "_conf_low").alias("conf_low"),
                pl.col(ep + "_conf_high").alias("conf_high"),
            )
        )
cols = ["compound_id", "endpoint", "value", "conf_low", "conf_high"]
a = pl.concat(raw_labels).unique().sort(cols)
b = (
    measurements.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))
    .select(cols)
    .unique()
    .sort(cols)
)
assert a.equals(b)
oof = read("tabicl-oof.parquet")
z = oof.join(b, on=["compound_id", "endpoint"], how="left", validate="1:1", suffix="_raw")
assert z["value"].null_count() == 0
for old, new in [("y_true", "value"), ("conf_low", "conf_low_raw"), ("conf_high", "conf_high_raw")]:
    np.testing.assert_allclose(z[old], z[new], rtol=0, atol=0, equal_nan=True)
summary["labels"] = dict(raw_direct_labels=a.height, oof_labels_checked=z.height, exact_match=True)

submission = read("regression.parquet")
source = read("tabicl-regression.parquet")
test = pl.read_csv(path("cyp-challenge-TEST-BLINDED.csv"))
assert (
    submission.select("Molecule_Name", "SMILES")
    .sort("Molecule_Name")
    .equals(test.select("Molecule_Name", "SMILES").sort("Molecule_Name"))
)
assert submission.sort("Molecule_Name").equals(
    source.select(submission.columns).sort("Molecule_Name")
)
folds = read("tabicl-test-folds.parquet")
means = (
    folds.group_by("compound_id", "endpoint")
    .agg(pl.col("prediction").mean())
    .pivot(on="endpoint", index="compound_id", values="prediction")
)
joined = source.join(means, on="compound_id", validate="1:1", suffix="_recomputed")
for ep in [c for c in source.columns if c.endswith("_pIC50_direct_inhibition")]:
    np.testing.assert_allclose(joined[ep], joined[ep + "_recomputed"], atol=1e-12, rtol=0)
summary["submission"] = dict(
    rows=submission.height,
    exact_raw_smiles=True,
    packaging_unchanged=True,
    fold_mean_reproduced=True,
)
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2), flush=True)
