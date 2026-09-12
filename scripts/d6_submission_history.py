"""Existing receipts only: check whether our scored submissions identify blind moments."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import pearsonr, spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


frames = [pl.read_parquet(path(i)).sort("Molecule_Name") for i in range(3)]
assert all(f.height == 750 and f["Molecule_Name"].n_unique() == 750 for f in frames)
for f in frames[1:]:
    for c in ["Molecule_Name", "SMILES"]:
        assert f[c].equals(frames[0][c])
rows = []
vecs = []
for i, item in enumerate(cfg["submissions"]):
    intent = json.loads(path(i + 6).read_text())
    assert intent["submission_id"] == item["submission_id"]
    assert hashlib.sha256(path(i).read_bytes()).hexdigest() == intent["file_sha256"]
    lb = pl.read_csv(path(i + 3))
    r = lb.filter((pl.col("Username") == "N283T") & (pl.col("Submitted") == item["submitted"]))
    assert r.height == 1
    r = r.row(0, named=True)
    official = lb.filter(pl.col("Username") == "TabICL-baseline")
    assert official.height == 1
    p = frames[i][cfg["endpoint"]].to_numpy().astype(float)
    assert np.isfinite(p).all()
    vecs.append(p)
    row = dict(
        name=item["name"],
        submission_id=item["submission_id"],
        submitted=item["submitted"],
        n=750,
        pred_mean=float(p.mean()),
        pred_sd=float(p.std(ddof=1)),
        pred_min=float(p.min()),
        pred_max=float(p.max()),
        pred_q10=float(np.quantile(p, 0.1)),
        pred_q90=float(np.quantile(p, 0.9)),
        LB_ST_RAE=float(r["ST-RAE"]),
        LB_MAE=float(r["MAE"]),
        LB_R2=float(r["R²"]),
        LB_Spearman=float(r["Spearman's ρ"]),
        LB_ST_RAE_std=float(r["ST-RAE std"]),
        LB_Spearman_std=float(r["Spearman's ρ std"]),
    )
    if official.height == 1:
        o = official.row(0, named=True)
        row.update(
            official_ST_RAE=float(o["ST-RAE"]),
            official_Spearman=float(o["Spearman's ρ"]),
            ST_RAE_gap=float(r["ST-RAE"] - o["ST-RAE"]),
            Spearman_gap=float(r["Spearman's ρ"] - o["Spearman's ρ"]),
        )
    rows.append(row)
pairs = []
for i in range(3):
    for j in range(i + 1, 3):
        x, y = vecs[i], vecs[j]
        slope, intercept = np.linalg.lstsq(np.column_stack([x, np.ones(len(x))]), y, rcond=None)[0]
        res = y - (slope * x + intercept)
        pairs.append(
            dict(
                source=cfg["submissions"][i]["name"],
                target=cfg["submissions"][j]["name"],
                slope=float(slope),
                intercept=float(intercept),
                Pearson=float(pearsonr(x, y).statistic),
                Spearman=float(spearmanr(x, y).statistic),
                affine_residual_RMSE=float(np.sqrt(np.mean(res**2))),
                affine_max_abs_residual=float(abs(res).max()),
                affine_exact_at_1e_8=bool(abs(res).max() < 1e-8),
            )
        )
pl.DataFrame(rows).write_csv(out / "submission-comparison.csv")
pl.DataFrame(pairs).write_csv(out / "affine-checks.csv")
frames[0].select("Molecule_Name", "SMILES").with_columns(
    [pl.Series(item["name"], p) for item, p in zip(cfg["submissions"], vecs, strict=True)]
).write_parquet(out / "aligned-test-predictions.parquet")
conclusion = dict(
    exact_receipt_hashes=True,
    exact_750_ids_smiles=True,
    exact_submission_timestamps=True,
    all_pairs_affine=all(p["affine_exact_at_1e_8"] for p in pairs),
    identifiable=False,
    reason=(
        "3 different prediction vectors require 3 independent label-prediction covariances "
        "in addition to unknown label mean and variance: 3 R2 equations vs at least 5 unknowns, "
        "even if scored-half prediction moments were known. Actual scored-half membership "
        "is unavailable. Spearman is not Pearson. Full750 moments differ from scored-half moments."
    ),
    no_network=True,
    no_inferred_blind_moments=True,
)
(out / "identifiability.json").write_text(json.dumps(conclusion, indent=2))
print(pl.DataFrame(rows).write_csv())
print(pl.DataFrame(pairs).write_csv())
print(json.dumps(conclusion))
