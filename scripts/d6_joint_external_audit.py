"""Independent identity, split, scaler, gradient and metric audit of joint graph fits."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


def read(i):
    return pl.read_parquet(path(i))


if cfg.get("outer_refit"):
    selected = {(r["fold"], r["method"]): r["best_epoch"] for r in json.loads(path(25).read_text())}
    for r in json.loads(path(24).read_text()):
        assert r["epochs"] == r["best_epoch"] == selected[(r["fold"], r["method"])]
        assert r["inner_ST_RAE"] is None


oof, members, external_members, grad, labels, assignments, compounds, ext = [
    pl.read_csv(path(i)) if i == 3 else read(i) for i in range(8)
]
keys = ["compound_id", "endpoint", "fold"]
ep = "CYP2D6_pIC50_direct_inhibition"
expected = (
    labels.filter(pl.col("endpoint") == ep)
    .join(assignments, on="compound_id", validate="1:1")
    .filter(~pl.col("is_anchor"))
    .sort(keys)
)
assert expected.height == 1457
spec = importlib.util.spec_from_file_location("official", path(8))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
forbidden = {
    Chem.MolToInchiKey(Chem.MolFromSmiles(s)).split("-")[0]
    for s in compounds["standardized_smiles"]
}
assert set(external_members["connectivity"]).isdisjoint(forbidden)
assert set(external_members["compound_id"]) == set(ext["compound_id"])
external_y = ext.select("Max_Response", "Activity at 57.14 uM", "inactive").to_numpy().astype(float)
external_y[:, :2] = np.clip(external_y[:, :2], -150, 50)
labeltable = labels.filter(pl.col("endpoint") == ep).select("compound_id", "value")
checks = []
for fold in range(5):
    m = members.filter(pl.col("outer_fold") == fold).join(
        assignments, on="compound_id", how="left", validate="1:1"
    )
    outer = (m["fold"].fill_null(-1).to_numpy() == fold) & ~m["is_anchor"].fill_null(
        False
    ).to_numpy()
    assert np.array_equal(outer, m["outer_holdout"].to_numpy())
    tr = m.filter("training")
    vi = m.filter("inner_validation")
    assert not set(tr["compound_id"]) & set(vi["compound_id"])
    assert not set(tr["compound_id"]) & set(m.filter("outer_holdout")["compound_id"])
    assert set(tr["group"]).isdisjoint(set(vi["group"]))
    yy = (
        tr.select("compound_id")
        .join(labeltable, on="compound_id", validate="1:1")["value"]
        .to_numpy()
    )
    if cfg.get("outer_refit"):
        required = labeltable.join(assignments, on="compound_id", validate="1:1").filter(
            (pl.col("fold") != fold) | pl.col("is_anchor")
        )
        assert set(tr["compound_id"]) == set(required["compound_id"])
        assert vi.is_empty()
    means = np.array([yy.mean(), *np.nanmean(external_y[:, :2], 0), 0.0])
    stds = np.array([yy.std(), *np.nanstd(external_y[:, :2], 0), 1.0])
    stds = np.maximum(stds, 0.05)
    counts = np.array([len(yy), *np.isfinite(external_y).sum(0)])
    for j, weight in enumerate([0, 0.05, 0.30]):
        s = json.loads(path(9 + fold * 3 + j).read_text())
        assert np.allclose(s["means"], means, rtol=1e-5, atol=1e-5)
        assert np.allclose(s["stds"], stds, rtol=1e-5, atol=1e-5)
        assert np.array_equal(s["observed_label_counts"], counts)
        assert np.allclose(s["weighted_label_mass"], counts * np.array([1, weight, weight, weight]))
        row = grad.filter((pl.col("fold") == fold) & (pl.col("public_weight") == weight)).row(
            0, named=True
        )
        assert row["direct_mp_gradient_norm"] > 0
        assert (
            row["public_mp_gradient_norm"] == 0
            if weight == 0
            else row["public_mp_gradient_norm"] > 0
        )
        checks.append(
            dict(
                fold=fold,
                public_weight=weight,
                training_n=len(yy),
                scaler_verified=True,
                encoder_public_gradient_verified=True,
            )
        )
rows = []
for (method,), g in oof.group_by("method"):
    g = g.sort(keys)
    assert g.height == 1457 and g.select(keys).n_unique() == 1457
    for col in keys + ["conf_low", "conf_high", "is_neighbor"]:
        assert g[col].equals(expected[col])
    assert g["y_true"].equals(expected["value"])
    assert g["prediction"].is_finite().all()
    for subset in ["all", "neighbors"]:
        h = g if subset == "all" else g.filter("is_neighbor")
        rows.append(
            dict(
                method=method,
                subset=subset,
                ST_RAE=float(
                    scorer.rae_soft_threshold_absolute_error(
                        h["y_true"].to_numpy(),
                        h["prediction"].to_numpy(),
                        y_true_lower=h["conf_low"].to_numpy(),
                        y_true_upper=h["conf_high"].to_numpy(),
                    )
                ),
            )
        )
pl.DataFrame(rows).write_csv(out / "verified-metrics.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            checks=checks,
            methods=oof["method"].n_unique(),
            rows_per_method=1457,
            all_challenge_connectivity_excluded=True,
        ),
        indent=2,
    )
)
print(pl.DataFrame(rows))
