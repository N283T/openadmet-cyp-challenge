"""Independent saved-booster reconstruction and inner-MCC threshold audit."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.metrics import matthews_corrcoef

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
cfg = json.loads((inp / config["source_config"]).read_text())


def path(name):
    found = [p for p in inp.iterdir() if p.name.endswith("-" + name)]
    assert len(found) == 1, (name, found)
    return found[0]


def read(key):
    return pl.read_parquet(path(cfg[key]))


comp = read("compounds").sort("compound_id")
ids = comp.select("compound_id")
index = {c: i for i, c in enumerate(comp["compound_id"])}
labels = read("labels").filter(pl.col("endpoint").is_in(["CYP2D6_is_TDI", "CYP3A4_is_TDI"]))
assign = (
    read("assignments").filter(pl.col("strategy") == "umap").select("compound_id", "fold", "group")
)
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
assert labels.height == 5081 and test.height == 750
assert set(test["compound_id"]).isdisjoint(set(assign["compound_id"]))
oof = pl.read_parquet(path("oof.parquet"))
blind = pl.read_parquet(path("fold-test-probabilities.parquet"))
members = pl.read_parquet(path("membership.parquet"))
models = json.loads(path("models.json").read_text())
assert oof.height == 5081 * 12 and blind.height == 750 * 2 * 5 * 12 and len(models) == 120
assert not oof.select("compound_id", "endpoint", "method").is_duplicated().any()
for (endpoint, outer), m in members.group_by("endpoint", "outer_fold"):
    expected = (
        labels.filter(pl.col("endpoint") == endpoint)
        .join(assign, on="compound_id")
        .sort("compound_id")
    )
    m = m.sort("compound_id")
    assert m.select("compound_id", "fold").equals(expected.select("compound_id", "fold"))
    assert np.array_equal(m["outer_validation"].to_numpy(), m["fold"].to_numpy() == outer)
    assert np.array_equal(m["inner_validation"].to_numpy(), m["fold"].to_numpy() == (outer + 1) % 5)
    assert not m.filter(pl.col("outer_validation") & pl.col("inner_validation")).height


def matrix(table, prefix=None):
    columns = [
        c for c in table.columns if c != "compound_id" and (prefix is None or c.startswith(prefix))
    ]
    x = (
        ids.join(
            table.select("compound_id", *columns),
            on="compound_id",
            how="left",
            validate="1:1",
            maintain_order="left",
        )
        .select(columns)
        .to_numpy()
        .astype(np.float32)
    )
    x[~np.isfinite(x)] = np.nan
    return x


xm = matrix(pl.concat([pl.read_parquet(path(n)) for n in cfg["mordred"]]))
xc = matrix(read("class_emb"), "emb_")
xx = matrix(read("mixed_emb"), "emb_")
fingerprinter = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
xf = np.stack(
    [
        fingerprinter.GetFingerprintAsNumPy(Chem.MolFromSmiles(s))
        for s in comp["standardized_smiles"]
    ]
).astype(np.float32)
views = {
    "morgan": xf,
    "mordred": xm,
    "fp_mordred": np.column_stack([xf, xm]),
    "external_class": xc,
    "external_mixed": xx,
    "external_class_mordred": np.column_stack([xc, xm]),
}
inner_rows, checks = [], []
for record in models:
    ep, method, outer = [record[k] for k in ["endpoint", "method", "fold"]]
    view = method.rsplit("_", 1)[0]
    x = views[view]
    h = labels.filter(pl.col("endpoint") == ep).join(assign, on="compound_id").sort("compound_id")
    ff = h["fold"].to_numpy()
    inner = (outer + 1) % 5
    assert record["n_train"] == int((ff != outer).sum())
    assert record["n_inner_fit"] == int(((ff != outer) & (ff != inner)).sum())
    assert record["n_inner_validation"] == int((ff == inner).sum())
    assert record["dimensions"] == x.shape[1]
    hold = h.filter(pl.col("fold") == inner)
    ix = np.array([index[c] for c in hold["compound_id"]])
    booster = lgb.Booster(model_file=str(path(f"{method}-{ep}-{outer}-inner.txt")))
    p = booster.predict(x[ix], num_threads=2)
    y = hold["value"].to_numpy().astype(int)
    # Independent exhaustive threshold contingency matrices, rather than source cumulative sort.
    candidates = np.unique(np.r_[p, 0.5])
    calls = p[:, None] >= candidates[None, :]
    tp = (calls * y[:, None]).sum(axis=0)
    fp = calls.sum(axis=0) - tp
    fn = y.sum() - tp
    tn = len(y) - tp - fp - fn
    den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    scores = np.divide(tp * tn - fp * fn, den, out=np.zeros_like(den, dtype=float), where=den != 0)
    selected = max(zip(scores, -np.abs(candidates - 0.5), candidates, strict=True))[2]
    assert abs(selected - record["threshold"]) < 1e-12
    inner_rows.append(
        hold.select(
            "compound_id",
            "endpoint",
            pl.col("fold").alias("inner_fold"),
            "group",
            pl.col("value").alias("y_true"),
        ).with_columns(
            pl.lit(outer).alias("outer_fold"),
            pl.lit(method).alias("method"),
            pl.Series("prediction", p),
            pl.lit(float(selected)).alias("threshold"),
            pl.lit(cfg["split_run"]).alias("split_id"),
        )
    )
    outermodel = lgb.Booster(model_file=str(path(f"{method}-{ep}-{outer}.txt")))
    errors = []
    for is_test, table in [(False, oof), (True, blind)]:
        saved = table.filter(
            (pl.col("method") == method) & (pl.col("endpoint") == ep) & (pl.col("fold") == outer)
        ).sort("compound_id")
        expected = test if is_test else h.filter(pl.col("fold") == outer).select("compound_id")
        assert saved["compound_id"].equals(expected.sort("compound_id")["compound_id"])
        assert (saved["split_id"] == cfg["split_run"]).all()
        if not is_test:
            assert np.array_equal(
                saved["y_true"].to_numpy(), h.filter(pl.col("fold") == outer)["value"].to_numpy()
            )
        recovered = outermodel.predict(
            x[np.array([index[c] for c in saved["compound_id"]])], num_threads=2
        )
        error = float(np.max(np.abs(recovered - saved["prediction"].to_numpy())))
        assert error < 1e-12 and (saved["threshold"] == record["threshold"]).all()
        errors.append(error)
    checks.append(
        dict(
            endpoint=ep,
            method=method,
            fold=outer,
            inner_n=len(y),
            threshold=float(selected),
            inner_MCC=float(matthews_corrcoef(y, p >= selected)),
            outer_error=errors[0],
            test_error=errors[1],
        )
    )
    print(ep, method, outer, flush=True)
inner = pl.concat(inner_rows)
assert inner.height == 5081 * 12
inner.write_parquet(out / "inner-predictions.parquet")
pl.DataFrame(checks).write_csv(out / "reconstruction.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            status="passed",
            oof_rows=oof.height,
            test_fold_rows=blind.height,
            inner_rows=inner.height,
            models=len(checks),
            thresholds="All120 independently reconstructed exhaustive MCC optima",
            membership="Original UMAP exact; provisional training excludes outer and inner folds",
            prediction="AllouterOOFand750test/fold/endpoint/method match savedboosterswithin1e-12",
        ),
        indent=2,
    )
)
