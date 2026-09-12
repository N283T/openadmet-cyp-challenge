"""Prepare exact old-A4 restoration and a fixed single-LGBM D6 alternative."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from cyp_challenge.tracks import ENDPOINTS, validate

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(k):
    return pl.read_parquet(inp / cfg[k])


reference = pl.read_csv(inp / cfg["reference"])
aliases = (
    read("aliases")
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
identities = reference.select("SMILES", "Molecule_Name").join(
    aliases, on="Molecule_Name", validate="1:1"
)
assert identities.height == identities["compound_id"].n_unique() == 750
old = read("old_submission").select("SMILES", "Molecule_Name", "CYP3A4_is_TDI")
assert (
    old.select("SMILES", "Molecule_Name")
    .sort("Molecule_Name")
    .equals(reference.select("SMILES", "Molecule_Name").sort("Molecule_Name"))
)
mean = read("current_test").filter(
    (pl.col("endpoint") == "CYP2D6_is_TDI") & (pl.col("method") == "all_plain_pfn_graph_mean")
)
single = read("lgbm_test").filter(
    (pl.col("endpoint") == "CYP2D6_is_TDI") & (pl.col("method") == "morgan_balanced")
)
assert single.height == 3750
single = single.group_by("compound_id").agg(
    pl.col("prediction").mean(), pl.col("fold").n_unique().alias("folds")
)
assert (single["folds"] == 5).all()
single = single.with_columns(
    pl.lit(0.5).alias("threshold"), (pl.col("prediction") >= 0.5).alias("binary")
)
reports = []
for name, t in [("retain_d6_restore_a4", mean), ("compact_morgan_restore_a4", single)]:
    wide = identities.join(
        t.select("compound_id", pl.col("binary").alias("CYP2D6_is_TDI")),
        on="compound_id",
        validate="1:1",
    ).join(old, on=["SMILES", "Molecule_Name"], validate="1:1")
    wide = wide.select("SMILES", "Molecule_Name", *ENDPOINTS["classification"]).sort(
        "Molecule_Name"
    )
    assert wide.height == 750
    assert wide["CYP3A4_is_TDI"].equals(old.sort("Molecule_Name")["CYP3A4_is_TDI"])
    wide.write_parquet(out / f"{name}.parquet")
    wide.write_csv(out / f"{name}.csv")
    validation = validate(out / f"{name}.parquet", "classification", inp / cfg["reference"])
    reports.append(
        dict(
            name=name,
            validation=validation,
            positives={ep: int(wide[ep].sum()) for ep in ENDPOINTS["classification"]},
        )
    )
a = (
    read("lgbm_oof")
    .filter((pl.col("endpoint") == "CYP2D6_is_TDI") & (pl.col("method") == "morgan_balanced"))
    .sort("compound_id")
)
b = (
    read("current_oof")
    .filter(
        (pl.col("endpoint") == "CYP2D6_is_TDI") & (pl.col("method") == "all_plain_pfn_graph_mean")
    )
    .sort("compound_id")
)
assert a["compound_id"].equals(b["compound_id"]) and np.array_equal(a["y_true"], b["y_true"])
y = a["y_true"].to_numpy()
pa, pb = (a["prediction"] >= 0.5).to_numpy(), (b["prediction"] >= b["threshold"]).to_numpy()
groups = a["group"].to_numpy()
indices = [np.flatnonzero(groups == g) for g in np.unique(groups)]
rng, differences = np.random.default_rng(731), []
for _ in range(500):
    ix = np.concatenate([indices[i] for i in rng.integers(0, len(indices), len(indices))])
    differences.append(matthews_corrcoef(y[ix], pa[ix]) - matthews_corrcoef(y[ix], pb[ix]))
comparison = dict(
    single_MCC=float(matthews_corrcoef(y, pa)),
    reference_MCC=float(matthews_corrcoef(y, pb)),
    single_AUC=float(roc_auc_score(y, a["prediction"])),
    conditional_difference_lo=float(np.quantile(differences, 0.025)),
    conditional_difference_hi=float(np.quantile(differences, 0.975)),
    selection="Prior-screen fixed model; conditional comparison, not nested winner selection",
)
(out / "summary.json").write_text(
    json.dumps(dict(config=cfg, candidates=reports, comparison=comparison), indent=2)
)
print(json.dumps(dict(candidates=reports, comparison=comparison)))
