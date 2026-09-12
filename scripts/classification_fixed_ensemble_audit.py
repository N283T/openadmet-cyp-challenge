"""Independent equal-weight means, inner thresholds, metrics and deployment audit."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
audit = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
cfg = json.loads((inp / audit["config"]).read_text())


def path(name):
    files = [p for p in inp.iterdir() if p.name.endswith("-" + name)]
    assert len(files) == 1, (name, files)
    return files[0]


keys = ["compound_id", "endpoint", "fold"]
canonical = pl.read_parquet(inp / audit["canonical"])
saved = {k: pl.read_parquet(path("ensemble-" + k + ".parquet")) for k in ["oof", "inner", "test"]}
metrics = pl.read_csv(path("metrics.csv"))
thresholds = pl.read_csv(path("thresholds.csv"))
assert len(cfg["ensembles"]) == 7
assert saved["oof"].height == saved["inner"].height == 35567 and saved["test"].height == 10500
assert metrics.height == 168 and thresholds.height == 70
banks = {}
for source in cfg["sources"]:
    for kind in ["oof", "inner", "test"]:
        frame = pl.read_parquet(path(source[kind]))
        if kind == "inner":
            if "fold" in frame.columns:
                frame = frame.drop("fold")
            frame = frame.rename({"outer_fold": "fold"})
        prediction = "probability" if "probability" in frame.columns else "prediction"
        for (method,), t in frame.group_by("method"):
            label = source["prefix"] + "::" + method
            cols = keys + [pl.col(prediction).cast(pl.Float64).alias("prediction")]
            if kind != "test":
                cols += [pl.col("y_true").cast(pl.Int64)]
            t = t.select(cols).sort(keys)
            assert t.height == (7500 if kind == "test" else 5081)
            assert not t.select(keys).is_duplicated().any()
            if kind != "test":
                check = t.join(
                    canonical.select(
                        "compound_id",
                        "endpoint",
                        pl.col("fold").alias("canonical_fold"),
                        pl.col("y_true").cast(pl.Int64).alias("canonical_y"),
                    ),
                    on=["compound_id", "endpoint"],
                    validate="1:1",
                )
                assert check.height == 5081 and (check["y_true"] == check["canonical_y"]).all()
                expected = (check["fold"] + 1) % 5 if kind == "inner" else check["fold"]
                assert (check["canonical_fold"] == expected).all()
            banks[kind, label] = t


def threshold(y, p):
    cuts = np.unique(np.r_[p, 0.5])
    call = p[:, None] >= cuts[None, :]
    tp = (call * y[:, None]).sum(0)
    fp = call.sum(0) - tp
    fn = y.sum() - tp
    tn = len(y) - tp - fp - fn
    den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    score = np.divide(tp * tn - fp * fn, den, out=np.zeros_like(den, dtype=float), where=den > 0)
    return max(zip(score, -np.abs(cuts - 0.5), cuts, strict=True))[2]


checks = []
for name, members in cfg["ensembles"].items():
    reconstructed = {}
    for kind in ["oof", "inner", "test"]:
        parts = [banks[kind, member] for member in members]
        ref = parts[0]
        assert all(part.select(keys).equals(ref.select(keys)) for part in parts)
        # Independently aggregate a long table rather than source numpy list-average.
        combined = (
            pl.concat(parts, how="vertical_relaxed")
            .group_by(keys)
            .agg(pl.col("prediction").mean())
            .sort(keys)
        )
        if kind != "test":
            combined = combined.join(ref.select(keys + ["y_true"]), on=keys, validate="1:1")
            actual = saved[kind].filter(pl.col("method") == name).sort(keys)
            assert actual.select(keys).equals(combined.select(keys))
            error = float(
                np.max(np.abs(actual["prediction"].to_numpy() - combined["prediction"].to_numpy()))
            )
            assert error < 1e-12
        reconstructed[kind] = combined
    for ep in canonical["endpoint"].unique():
        cuts = []
        for fold in range(5):
            inner = reconstructed["inner"].filter(
                (pl.col("endpoint") == ep) & (pl.col("fold") == fold)
            )
            outer = reconstructed["oof"].filter(
                (pl.col("endpoint") == ep) & (pl.col("fold") == fold)
            )
            assert set(inner["compound_id"]).isdisjoint(set(outer["compound_id"]))
            cut = threshold(inner["y_true"].to_numpy(), inner["prediction"].to_numpy())
            recorded = thresholds.filter(
                (pl.col("method") == name) & (pl.col("endpoint") == ep) & (pl.col("fold") == fold)
            )
            # Floating arithmetic order may change an equivalent threshold by <1e-12.
            assert recorded.height == 1 and abs(cut - recorded["threshold"][0]) < 1e-12
            cuts.append(recorded["threshold"][0])
        deployed = (
            reconstructed["test"]
            .filter(pl.col("endpoint") == ep)
            .group_by("compound_id", "endpoint")
            .agg(pl.col("prediction").mean())
            .sort("compound_id")
        )
        actual = (
            saved["test"]
            .filter((pl.col("method") == name) & (pl.col("endpoint") == ep))
            .sort("compound_id")
        )
        assert actual.height == 750 and actual["compound_id"].equals(deployed["compound_id"])
        error = float(
            np.max(np.abs(deployed["prediction"].to_numpy() - actual["prediction"].to_numpy()))
        )
        assert error < 1e-12 and np.allclose(
            actual["threshold"], np.median(cuts), rtol=0, atol=1e-12
        )
        assert np.array_equal(actual["binary"], actual["prediction"] >= actual["threshold"])
        for fold in [-1, 0, 1, 2, 3, 4]:
            t = saved["oof"].filter((pl.col("method") == name) & (pl.col("endpoint") == ep))
            if fold != -1:
                t = t.filter(pl.col("fold") == fold)
            y, p = t["y_true"].to_numpy(), t["prediction"].to_numpy()
            for policy, cut in [("fixed_05", 0.5), ("inner_mcc", t["threshold"].to_numpy())]:
                row = metrics.filter(
                    (pl.col("method") == name)
                    & (pl.col("endpoint") == ep)
                    & (pl.col("fold") == fold)
                    & (pl.col("policy") == policy)
                )
                assert row.height == 1
                values = dict(
                    MCC=matthews_corrcoef(y, p >= cut),
                    AUROC=roc_auc_score(y, p),
                    AP=average_precision_score(y, p),
                )
                assert all(abs(float(row[c][0]) - v) < 1e-12 for c, v in values.items())
        checks.append(
            dict(
                method=name,
                endpoint=ep,
                test_rows=750,
                test_mean_error=error,
                median_threshold=float(np.median(cuts)),
            )
        )
pl.DataFrame(checks).write_csv(out / "reconstruction.csv")
metrics.filter(pl.col("fold") == -1).write_csv(out / "pooled-metrics.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            status="passed",
            ensembles=7,
            oof_rows=35567,
            inner_rows=35567,
            test_rows=10500,
            inner_thresholds=70,
            metrics_checked=168,
            source_memberships="CanonicalUMAPandlabels match; innerfold=(outer+1)%5",
            deployment="Equalmembermeans thenfivefoldmean; medianinnerthreshold; nofull-labelrefit",
        ),
        indent=2,
    )
)
