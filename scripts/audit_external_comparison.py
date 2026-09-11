"""Audit paired five-fold regression predictions and estimate cluster-resampled deltas."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

from cyp_challenge.storage import write_json
from cyp_challenge.tracks import ENDPOINTS


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    labels = pl.read_parquet(inp / cfg["labels"]).filter(
        pl.col("endpoint").is_in(ENDPOINTS["regression"])
    )
    splits = pl.read_parquet(inp / cfg["splits"]).filter(pl.col("strategy") == "umap")
    test_keys = (
        pl.read_parquet(inp / cfg["test_reference"])
        .filter(pl.col("endpoint").is_in(ENDPOINTS["regression"]))
        .select("compound_id", "endpoint")
    )
    assert test_keys.height == 3000 == test_keys.n_unique()
    spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
    official = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(official)
    variants, scores, memberships = {}, [], {}
    for variant, files in cfg["variants"].items():
        frames, tests = [], []
        for item in files:
            frame = pl.read_parquet(inp / item["oof"])
            assert frame.select("compound_id", "endpoint").n_unique() == frame.height
            if "summary" in item:
                summary = json.loads((inp / item["summary"]).read_text())
                assert summary["config"]["variant"] == variant
                member = pl.read_parquet(inp / item["membership"])
                memberships[(variant, summary["config"]["fold"])] = member
                test = pl.read_parquet(inp / item["test"])
                assert test.select("compound_id", "endpoint").n_unique() == test.height == 3000
                assert test_keys.join(test, on=["compound_id", "endpoint"]).height == 3000
                assert (test["fold"] == item["fold"]).all()
                assert np.isfinite(test["prediction"].to_numpy()).all()
                assert set(member["compound_id"]).isdisjoint(test["compound_id"])
                tests.append(test)
            frames.append(frame)
        if tests:
            all_test = pl.concat(tests)
            all_test.write_parquet(out / f"{variant}-test-members.parquet")
            ensemble = all_test.group_by("compound_id", "endpoint", "split_id").agg(
                pl.col("prediction").mean(), pl.col("fold").n_unique().alias("n_folds")
            )
            assert ensemble.height == 3000 and (ensemble["n_folds"] == 5).all()
            ensemble.write_parquet(out / f"{variant}-test-ensemble.parquet")
        pred = pl.concat(frames).filter(pl.col("endpoint").is_in(ENDPOINTS["regression"]))
        assert pred.height == labels.height == pred.select("compound_id", "endpoint").n_unique()
        f = labels.join(pred, on=["compound_id", "endpoint"], validate="1:1")
        assert f.height == labels.height and (f["value"] == f["y_true"]).all()
        f = f.join(
            splits.select("compound_id", pl.col("fold").alias("expected_fold"), "group"),
            on="compound_id",
            validate="m:1",
        )
        assert (f["fold"] == f["expected_fold"]).all()
        assert f["group"].null_count() == 0 and np.isfinite(f["prediction"].to_numpy()).all()
        variants[variant] = f
        f.with_columns(pl.lit(variant).alias("variant")).write_parquet(
            out / f"{variant}-oof.parquet"
        )
        for ep in ENDPOINTS["regression"]:
            for fold in [-1, 0, 1, 2, 3, 4]:
                d = f.filter(pl.col("endpoint") == ep)
                if fold >= 0:
                    d = d.filter(pl.col("fold") == fold)
                y, p = d["value"].to_numpy(), d["prediction"].to_numpy()
                metrics = dict(
                    ST_RAE=float(
                        official.rae_soft_threshold_absolute_error(
                            y,
                            p,
                            y_true_upper=d["conf_high"].to_numpy(),
                            y_true_lower=d["conf_low"].to_numpy(),
                        )
                    ),
                    MAE=float(np.abs(p - y).mean()),
                    Spearman=float(spearmanr(y, p).statistic),
                )
                scores.append(dict(variant=variant, endpoint=ep, fold=fold, n=len(y), **metrics))
                if variant != "baseline" and fold >= 0:
                    summary = json.loads(
                        (inp / next(i["summary"] for i in files if i["fold"] == fold)).read_text()
                    )
                    saved = next(r for r in summary["metrics"] if r["endpoint"] == ep)
                    assert np.isclose(saved["ST-RAE"], metrics["ST_RAE"], atol=1e-10)
                    assert np.isclose(saved["MAE"], metrics["MAE"], atol=1e-10)
                    assert np.isclose(saved["Spearman"], metrics["Spearman"], atol=1e-10)
    for fold in range(5):
        a, b = memberships[("control", fold)], memberships[("external", fold)]
        assert a.sort("compound_id").equals(b.sort("compound_id"))
        assert set(a.filter(pl.col("role") == "train")["connectivity"]).isdisjoint(
            a.filter(pl.col("role").is_in(["inner_validation", "outer_holdout"]))["connectivity"]
        )
    pl.DataFrame(scores).write_csv(out / "metrics.csv")
    paired = variants["control"].join(
        variants["external"].select(
            "compound_id", "endpoint", pl.col("prediction").alias("external")
        ),
        on=["compound_id", "endpoint"],
        validate="1:1",
    )
    groups = sorted(paired["group"].unique())
    rng = np.random.default_rng(cfg["bootstrap_seed"])
    samples = rng.integers(0, len(groups), size=(cfg["bootstrap_replicates"], len(groups)))
    intervals, deltas = [], []
    for ep in ENDPOINTS["regression"]:
        f = paired.filter(pl.col("endpoint") == ep)
        y, p, e = f["value"].to_numpy(), f["prediction"].to_numpy(), f["external"].to_numpy()
        lo, hi = f["conf_low"].to_numpy(), f["conf_high"].to_numpy()
        index = {g: np.flatnonzero(f["group"].to_numpy() == g) for g in groups}
        delta = []
        for draw in samples:
            ix = np.concatenate([index[groups[i]] for i in draw])
            denominator = (
                np.maximum(y[ix].mean() - hi[ix], 0) + np.maximum(lo[ix] - y[ix].mean(), 0)
            ).sum()
            a = np.maximum(p[ix] - hi[ix], 0) + np.maximum(lo[ix] - p[ix], 0)
            b = np.maximum(e[ix] - hi[ix], 0) + np.maximum(lo[ix] - e[ix], 0)
            delta.append(float((b.sum() - a.sum()) / denominator))
        deltas.append(delta)
        intervals.append(
            dict(
                endpoint=ep,
                delta_ST_RAE_q025=float(np.quantile(delta, 0.025)),
                delta_ST_RAE_q975=float(np.quantile(delta, 0.975)),
            )
        )
    macro = np.mean(deltas, axis=0)
    intervals.append(
        dict(
            endpoint="macro",
            delta_ST_RAE_q025=float(np.quantile(macro, 0.025)),
            delta_ST_RAE_q975=float(np.quantile(macro, 0.975)),
        )
    )
    pl.DataFrame(intervals).write_csv(out / "paired-cluster-bootstrap.csv")
    pooled = pl.DataFrame(scores).filter(pl.col("fold") == -1)
    macro_scores = (
        pooled.group_by("variant")
        .agg(pl.col("ST_RAE").mean(), pl.col("MAE").mean(), pl.col("Spearman").mean())
        .sort("variant")
    )
    macro_scores.write_csv(out / "macro.csv")
    write_json(
        out / "audit.json",
        dict(
            passed=True,
            config=cfg,
            macro=macro_scores.to_dicts(),
            pooled=pooled.to_dicts(),
            bootstrap=intervals,
            limitation=(
                "Paired UMAP cluster resampling of fixed OOF predictions; "
                "not training-seed uncertainty. No test outcomes used."
            ),
        ),
    )
    print(macro_scores, flush=True)


if __name__ == "__main__":
    main()
