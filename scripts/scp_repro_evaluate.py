"""Evaluate native community-recipe predictions and its prescribed direct average."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import r2_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


reference = pl.read_parquet(path(0)).filter(pl.col("method") == "caruana_weight_crossfit")
assignments = pl.read_parquet(path(1))
spec = importlib.util.spec_from_file_location("official", path(2))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
keys = ["compound_id", "endpoint", "fold"]
canonical = keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]
reference = reference.sort(keys)
assert len(reference) == 6384
checked = reference.join(
    assignments.select("compound_id", pl.col("fold").alias("assigned_fold"), "is_anchor"),
    on="compound_id",
    validate="m:1",
)
assert not checked.filter((pl.col("fold") != pl.col("assigned_fold")) | pl.col("is_anchor")).height
models = {}
tests = {}
provenance = []
for member in cfg["members"]:
    name = member["name"]
    t = pl.read_parquet(path(member["oof_input"]))
    if "value" in t.columns and "y_true" not in t.columns:
        t = t.rename({"value": "y_true"})
    if member.get("oof_method"):
        t = t.filter(pl.col("method") == member["oof_method"])
    assert not t.select(keys).is_duplicated().any(), name
    aligned = []
    for endpoint in member["endpoints"]:
        h = t.filter(pl.col("endpoint") == endpoint).sort(keys)
        ref = reference.filter(pl.col("endpoint") == endpoint)
        assert len(h) == len(ref), (name, endpoint, len(h), len(ref))
        for c in keys:
            assert h[c].equals(ref[c]), (name, c)
        assert np.allclose(h["y_true"], ref["y_true"], rtol=0, atol=1e-9), name
        for c in ["conf_low", "conf_high"]:
            if c in h.columns:
                assert np.allclose(h[c], ref[c], rtol=0, atol=1e-9), (name, c)
        assert h["prediction"].is_finite().all()
        aligned.append(ref.select(canonical).with_columns(h["prediction"]))
    models[name] = pl.concat(aligned).with_columns(pl.lit(name).alias("method"))
    test = pl.read_parquet(path(member["test_input"]))
    if member.get("test_method"):
        test = test.filter(pl.col("method") == member["test_method"])
    # Input test predictions are already the source recipe's five-fold average.
    tk = ["compound_id", "endpoint"]
    assert not test.select(tk).is_duplicated().any(), name
    for endpoint in member["endpoints"]:
        assert test.filter(pl.col("endpoint") == endpoint).height == 750, (name, endpoint)
    assert test["prediction"].is_finite().all()
    tests[name] = test.select(tk + ["prediction"]).with_columns(pl.lit(name).alias("method"))
    provenance.append(member)
means, blindmeans = [], []
for endpoint in sorted(reference["endpoint"].unique()):
    names = [m["name"] for m in cfg["members"] if endpoint in m["endpoints"]]
    expected = 4 if "2D6" in endpoint else 2
    assert len(names) == expected, (endpoint, names)
    parts = [models[n].filter(pl.col("endpoint") == endpoint).sort(keys) for n in names]
    for p in parts[1:]:
        for c in canonical:
            assert p[c].equals(parts[0][c]), (endpoint, c)
    means.append(
        parts[0].with_columns(
            pl.Series("prediction", np.mean([p["prediction"].to_numpy() for p in parts], axis=0)),
            pl.lit("community_direct_mean").alias("method"),
        )
    )
    parts = [tests[n].filter(pl.col("endpoint") == endpoint).sort(tk) for n in names]
    for p in parts[1:]:
        assert p["compound_id"].equals(parts[0]["compound_id"])
    blindmeans.append(
        parts[0].with_columns(
            pl.Series("prediction", np.mean([p["prediction"].to_numpy() for p in parts], axis=0)),
            pl.lit("community_direct_mean").alias("method"),
        )
    )
mean = pl.concat(means).sort(keys)
blind = pl.concat(blindmeans).sort(tk)
mean.write_parquet(out / "community-ensemble-oof.parquet")
blind.write_parquet(out / "community-ensemble-test.parquet")
all_oof = pl.concat(
    [
        reference.select(canonical + ["prediction"]).with_columns(
            pl.lit("previous_global").alias("method")
        ),
        *models.values(),
        mean,
    ],
    how="diagonal_relaxed",
)
all_oof.write_parquet(out / "all-model-oof.parquet")
metrics = []
for (method, endpoint), t in all_oof.group_by("method", "endpoint"):
    for fold in [-1, 0, 1, 2, 3, 4]:
        h = t if fold == -1 else t.filter(pl.col("fold") == fold)
        y, p = [h[c].to_numpy() for c in ["y_true", "prediction"]]
        metrics.append(
            dict(
                method=method,
                endpoint=endpoint,
                fold=fold,
                n=len(h),
                ST_RAE=float(
                    scorer.rae_soft_threshold_absolute_error(
                        y,
                        p,
                        y_true_lower=h["conf_low"].to_numpy(),
                        y_true_upper=h["conf_high"].to_numpy(),
                    )
                ),
                ST_RAE_author_local=float(
                    (
                        np.maximum(h["conf_low"].to_numpy() - p, 0)
                        + np.maximum(p - h["conf_high"].to_numpy(), 0)
                    ).sum()
                    / np.abs(y - y.mean()).sum()
                ),
                MAE=float(np.mean(np.abs(y - p))),
                Spearman=float(spearmanr(y, p).statistic),
                Pearson=float(pearsonr(y, p).statistic),
                R2=float(r2_score(y, p)),
                prediction_mean=float(p.mean()),
                prediction_std=float(p.std()),
                target_mean=float(y.mean()),
                target_std=float(y.std()),
            )
        )
metrics = pl.DataFrame(metrics).sort("method", "endpoint", "fold")
metrics.write_csv(out / "metrics.csv")
metrics.filter(pl.col("fold") == -1).group_by("method").agg(
    pl.col("ST_RAE").mean(),
    pl.col("ST_RAE_author_local").mean(),
    pl.col("endpoint").n_unique().alias("endpoints"),
).write_csv(out / "macro.csv")
pl.concat([*tests.values(), blind]).group_by("method", "endpoint").agg(
    pl.len().alias("n"),
    pl.col("prediction").mean().alias("mean"),
    pl.col("prediction").std(ddof=0).alias("std"),
    pl.col("prediction").min().alias("min"),
    pl.col("prediction").max().alias("max"),
).write_csv(out / "test-distributions.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            members=provenance,
            canonical_oof=6384,
            canonical_test=3000,
            recipe="equal direct predictions; four members for2D6, two forotherCYPs",
            validation=cfg["validation_description"],
            calibration="none",
        ),
        indent=2,
    )
)
print(metrics.filter(pl.col("fold") == -1).write_csv(), flush=True)
