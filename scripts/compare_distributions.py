"""Compare fixed distribution summaries and holdout-fitted corrections/quantile choices."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from calibration_screen import apply_map, fit_map, soft
from rdkit import Chem
from scipy.stats import spearmanr

from cyp_challenge.calibration_alignment import aligned_matrix
from cyp_challenge.external_training import snapshot_input
from cyp_challenge.storage import write_json


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def table(name):
        return pl.read_parquet(inp / name)

    labels = pl.read_parquet(snapshot_input(inp, "regression-oof-diagnostics.parquet"))
    baseline_oof = pl.concat([table(v["oof"]) for v in cfg["baseline"]]).with_columns(
        pl.lit("baseline").alias("member")
    )
    baseline_test = pl.concat([table(v["test"]) for v in cfg["baseline"]]).with_columns(
        pl.lit("baseline").alias("member")
    )
    frames = [
        baseline_oof.select("compound_id", "endpoint", "fold", "estimator", "prediction", "member")
    ]
    tests = [
        baseline_test.select("compound_id", "endpoint", "fold", "estimator", "prediction", "member")
    ]
    for source in cfg.get("additional", []):
        frames.append(table(source["oof"]).select(frames[0].columns))
        tests.append(table(source["test"]).select(tests[0].columns))
    oof = pl.concat(frames)
    test = pl.concat(tests)
    assert (
        oof.select("compound_id", "endpoint", "fold", "estimator", "member").n_unique()
        == oof.height
    )
    assert (
        test.select("compound_id", "endpoint", "fold", "estimator", "member").n_unique()
        == test.height
    )
    test = test.group_by("compound_id", "endpoint", "estimator", "member").agg(
        pl.col("prediction").mean(), pl.col("fold").n_unique().alias("nfold")
    )
    assert (test["nfold"] == 5).all()
    test = test.drop("nfold")
    estimators = cfg["estimators"]
    groups = {"baseline": ["baseline"]}
    if cfg.get("additional"):
        groups.update(
            compact=["compact"],
            pretrained=["pretrained"],
            blend3=["baseline", "compact", "pretrained"],
        )
    splits = pl.read_parquet(snapshot_input(inp, "assignments.parquet")).filter(
        pl.col("strategy") == "umap"
    )
    groupmap = dict(splits.select("compound_id", "group").iter_rows())
    compounds = pl.read_parquet(snapshot_input(inp, "compounds.parquet"))
    keymap = {
        cid: Chem.MolToInchiKey(Chem.MolFromSmiles(s)).split("-")[0]
        for cid, s in compounds.select("compound_id", "standardized_smiles").iter_rows()
    }
    spec = importlib.util.spec_from_file_location(
        "official", snapshot_input(inp, "custom_scoring_functions.py")
    )
    official = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(official)
    scores, outputs, testoutputs, params, selected = [], [], [], [], []
    bootdata = {}
    for ep in sorted(labels["endpoint"].unique()):
        f = labels.filter(pl.col("endpoint") == ep).sort("compound_id")
        ids = f["compound_id"].to_numpy()
        y, lo, hi, base = [
            f[c].to_numpy().astype(float) for c in ["y_true", "conf_low", "conf_high", "prediction"]
        ]
        folds = f["fold"].to_numpy()
        clusters = np.array([groupmap[c] for c in ids])
        keys = np.array([keymap[c] for c in ids])
        bootdata[ep] = dict(y=y, lo=lo, hi=hi, base=base, clusters=clusters, pred={})

        def record(name, q, repeat=-1, y=y, lo=lo, hi=hi, folds=folds, ep=ep, f=f):
            assert np.isfinite(q).all()
            for fold in [-1, 0, 1, 2, 3, 4]:
                ix = np.ones(len(y), dtype=bool) if fold == -1 else folds == fold
                score = float(
                    official.rae_soft_threshold_absolute_error(
                        y[ix], q[ix], y_true_lower=lo[ix], y_true_upper=hi[ix]
                    )
                )
                scores.append(
                    dict(
                        endpoint=ep,
                        method=name,
                        repeat=repeat,
                        fold=fold,
                        ST_RAE=score,
                        MAE=float(np.abs(y[ix] - q[ix]).mean()),
                        Spearman=float(spearmanr(y[ix], q[ix]).statistic),
                    )
                )
            outputs.append(
                f.select("compound_id", "endpoint", "fold", "split_id", "y_true").with_columns(
                    pl.Series("prediction", q),
                    pl.lit(name).alias("method"),
                    pl.lit(repeat).alias("repeat"),
                )
            )

        for group, members in groups.items():
            ff = oof.filter(
                (pl.col("endpoint") == ep)
                & pl.col("member").is_in(members)
                & pl.col("estimator").is_in(estimators)
            )
            ff = ff.group_by("compound_id", "fold", "estimator").agg(
                pl.col("prediction").mean(), pl.col("member").n_unique().alias("n")
            )
            assert (ff["n"] == len(members)).all()
            ff = ff.pivot(on="estimator", index=["compound_id", "fold"], values="prediction")
            ff = (
                f.select("compound_id", "fold")
                .join(ff, on=["compound_id", "fold"], validate="1:1")
                .sort("compound_id")
            )
            assert ff["compound_id"].to_list() == f["compound_id"].to_list()
            mat = aligned_matrix(f, ff, estimators)
            assert np.isfinite(mat).all()
            if group == "baseline":
                np.testing.assert_allclose(
                    mat[:, estimators.index("mean")], base, rtol=1e-5, atol=1e-5
                )
            tt = (
                test.filter(
                    (pl.col("endpoint") == ep)
                    & pl.col("member").is_in(members)
                    & pl.col("estimator").is_in(estimators)
                )
                .group_by("compound_id", "estimator")
                .agg(pl.col("prediction").mean())
                .pivot(on="estimator", index="compound_id", values="prediction")
                .sort("compound_id")
            )
            tm = tt.select(estimators).to_numpy().astype(float)
            assert tt.height == 750 and np.isfinite(tm).all()
            fixed = {k: mat[:, i] for i, k in enumerate(estimators)}
            fixedtest = {k: tm[:, i] for i, k in enumerate(estimators)}
            fixed["mean_median50"] = (fixed["mean"] + fixed["median"]) / 2
            fixedtest["mean_median50"] = (fixedtest["mean"] + fixedtest["median"]) / 2
            for estimator, q in fixed.items():
                name = group + ":" + estimator
                record(name, q)
                bootdata[ep]["pred"][name] = q
                testoutputs.append(
                    tt.select("compound_id").with_columns(
                        pl.lit(ep).alias("endpoint"),
                        pl.lit(name).alias("method"),
                        pl.lit(-1).alias("fold"),
                        pl.lit("test:mean-fold-point-estimates").alias("split_id"),
                        pl.Series("prediction", fixedtest[estimator]),
                    )
                )
            for repeat in range(cfg["repeats"]):
                qmaps = {
                    k: np.full(len(y), np.nan)
                    for k in [
                        "mean_st_half",
                        "median_st_half",
                        "mean_st_full",
                        "median_st_full",
                        "select_quantile",
                        "select_quantile_half",
                    ]
                }
                for fold in range(5):
                    available = np.unique(clusters[folds == fold])
                    shuffled = np.random.default_rng(cfg["seed"] + repeat * 100 + fold).permutation(
                        available
                    )
                    side = np.isin(clusters, shuffled[: len(shuffled) // 2])
                    for s in [False, True]:
                        va = (folds == fold) & (side == s)
                        tr = (folds == fold) & (side != s) & ~np.isin(keys, keys[va])
                        assert tr.sum() > 20 and va.sum() > 20
                        assert set(keys[tr]).isdisjoint(keys[va])
                        for estimator in ["mean", "median"]:
                            p = fixed[estimator]
                            m = fit_map(
                                "affine_st", p[tr], y[tr], lo[tr], hi[tr], np.ones(tr.sum())
                            )
                            qmaps[estimator + "_st_full"][va] = apply_map(m, p[va])
                            qmaps[estimator + "_st_half"][va] = (p[va] + apply_map(m, p[va])) / 2
                            params.append(
                                dict(
                                    endpoint=ep,
                                    group=group,
                                    estimator=estimator,
                                    repeat=repeat,
                                    fold=fold,
                                    half=int(s),
                                    model=m,
                                )
                            )
                        # Select prespecified summaries without evaluation labels.
                        loss = np.array(
                            [soft(mat[tr, i], lo[tr], hi[tr]).sum() for i in range(len(estimators))]
                        )
                        index = int(np.argmin(loss))
                        qmaps["select_quantile"][va] = mat[va, index]
                        qmaps["select_quantile_half"][va] = (mat[va, index] + fixed["mean"][va]) / 2
                        selected.append(
                            dict(
                                endpoint=ep,
                                group=group,
                                repeat=repeat,
                                fold=fold,
                                half=int(s),
                                estimator=estimators[index],
                            )
                        )
                for method, q in qmaps.items():
                    record(group + ":" + method, q, repeat)
            # Full-data maps produce test artifacts, without labeling their training score CV.
            for estimator in ["mean", "median"]:
                p = fixed[estimator]
                m = fit_map("affine_st", p, y, lo, hi, np.ones(len(p)))
                params.append(
                    dict(
                        endpoint=ep,
                        group=group,
                        estimator=estimator,
                        repeat=-1,
                        fold=-1,
                        half=-1,
                        model=m,
                    )
                )
                for strength, suffix in [(0.5, "half"), (1.0, "full")]:
                    tp = fixedtest[estimator]
                    q = (1 - strength) * tp + strength * apply_map(m, tp)
                    testoutputs.append(
                        tt.select("compound_id").with_columns(
                            pl.lit(ep).alias("endpoint"),
                            pl.lit(group + ":" + estimator + "_st_" + suffix).alias("method"),
                            pl.lit(-1).alias("fold"),
                            pl.lit("test:full-OOF-calibration").alias("split_id"),
                            pl.Series("prediction", q),
                        )
                    )
    score = pl.DataFrame(scores)
    macro = (
        score.filter(pl.col("fold") == -1)
        .group_by("method", "repeat")
        .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
    )
    summary = (
        macro.group_by("method")
        .agg(
            pl.col("ST_RAE", "MAE", "Spearman").mean(),
            pl.col("ST_RAE").min().alias("min"),
            pl.col("ST_RAE").max().alias("max"),
            pl.col("repeat").min().alias("repeat_min"),
        )
        .sort("ST_RAE")
    )
    # Fixed estimators only: shared cluster-bootstrap deltas vs original baseline mean.
    allclusters = sorted(set(g for d in bootdata.values() for g in d["clusters"]))
    draws = np.random.default_rng(cfg["seed"]).multinomial(
        len(allclusters), np.ones(len(allclusters)) / len(allclusters), size=1000
    )
    boots = {}
    for ep, d in bootdata.items():
        ww = draws[:, [allclusters.index(g) for g in d["clusters"]]]
        mu = ww @ d["y"] / ww.sum(axis=1)
        den = (
            ww
            * (
                np.maximum(d["lo"][None, :] - mu[:, None], 0)
                + np.maximum(mu[:, None] - d["hi"][None, :], 0)
            )
        ).sum(axis=1)
        for method, q in d["pred"].items():
            boots[ep, method] = (
                ww @ (soft(q, d["lo"], d["hi"]) - soft(d["base"], d["lo"], d["hi"])) / den
            )
    bootrows = []
    eps = list(bootdata)
    for method in bootdata[eps[0]]["pred"]:
        for ep in eps + ["macro"]:
            values = (
                np.mean([boots[e, method] for e in eps], axis=0)
                if ep == "macro"
                else boots[ep, method]
            )
            low, high = np.quantile(values, [0.025, 0.975])
            bootrows.append(
                dict(endpoint=ep, method=method, delta_low=float(low), delta_high=float(high))
            )
    score.write_csv(out / "distribution-metrics.csv")
    summary.write_csv(out / "distribution-summary.csv")
    macro.write_csv(out / "distribution-macro-repeats.csv")
    pl.DataFrame(bootrows).write_csv(out / "distribution-bootstrap.csv")
    pl.DataFrame(selected).write_csv(out / "quantile-selections.csv")
    pl.concat(outputs).write_parquet(out / "distribution-comparison-oof.parquet")
    pl.concat(testoutputs).write_parquet(out / "distribution-comparison-test.parquet")
    write_json(out / "parameters.json", params)
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            audit_passed=True,
            results=summary.to_dicts(),
            limitations=[
                "Fixed summaries use ordinary OOF; selected/calibrated summaries "
                "use small same-holdout cluster-half fitting.",
                "Test quantile output averages fold-level quantiles; "
                "these are not quantiles of the mixture distribution.",
                "Member/summary selection is exploratory. "
                "No hidden test score or selection-adjusted uncertainty.",
            ],
        ),
    )
    print(summary.head(20))


if __name__ == "__main__":
    main()
