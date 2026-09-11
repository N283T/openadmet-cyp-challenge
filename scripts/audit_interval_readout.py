"""Audit interval-target readouts and fixed blends against both existing baselines."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from calibration_screen import soft
from scipy.stats import spearmanr

from cyp_challenge.calibration_alignment import aligned_matrix
from cyp_challenge.external_training import snapshot_input
from cyp_challenge.storage import write_json


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    neighbor_oof = pl.read_parquet(snapshot_input(inp, "neighborhood-oof.parquet"))
    neighbor_test = pl.read_parquet(snapshot_input(inp, "neighborhood-test.parquet"))
    labels = pl.read_parquet(snapshot_input(inp, "regression-oof-diagnostics.parquet"))
    oo = pl.concat([pl.read_parquet(inp / r["oof"]) for r in cfg["runs"]])
    tt = pl.concat([pl.read_parquet(inp / r["test"]) for r in cfg["runs"]])
    test = tt.group_by("compound_id", "endpoint", "method").agg(
        pl.col("prediction").mean(), pl.col("fold").n_unique().alias("nfold")
    )
    assert (test["nfold"] == 5).all()
    reference = pl.read_parquet(snapshot_input(inp, "distribution-comparison-oof.parquet")).filter(
        (pl.col("repeat") == -1) & (pl.col("method") == "blend3:q0.55")
    )
    ref_test = pl.read_parquet(snapshot_input(inp, "distribution-comparison-test.parquet")).filter(
        pl.col("method").is_in(["baseline:mean", "blend3:q0.55"])
    )
    assignments = pl.read_parquet(snapshot_input(inp, "assignments.parquet")).filter(
        pl.col("strategy") == "umap"
    )
    groupmap = dict(assignments.select("compound_id", "group").iter_rows())
    spec = importlib.util.spec_from_file_location(
        "official", snapshot_input(inp, "custom_scoring_functions.py")
    )
    official = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(official)
    methods = sorted(oo["method"].unique())
    assert oo.select("compound_id", "endpoint", "method").n_unique() == oo.height
    assert oo.height == 6525 * len(methods)
    rows, preds, testpreds, bootdata = [], [], [], {}
    for ep in sorted(labels["endpoint"].unique()):
        f = labels.filter(pl.col("endpoint") == ep).sort("compound_id")
        y, lo, hi, p = [
            f[c].to_numpy().astype(float) for c in ["y_true", "conf_low", "conf_high", "prediction"]
        ]
        folds = f["fold"].to_numpy()
        proposed = reference.filter(pl.col("endpoint") == ep).rename({"prediction": "blend"})
        q55 = aligned_matrix(f, proposed, ["blend"])[:, 0]
        part = oo.filter(pl.col("endpoint") == ep).pivot(
            on="method", index=["compound_id", "fold"], values="prediction"
        )
        matrix = aligned_matrix(f, part, methods)
        t = (
            test.filter(pl.col("endpoint") == ep)
            .pivot(on="method", index="compound_id", values="prediction")
            .sort("compound_id")
        )
        rt = ref_test.filter(pl.col("endpoint") == ep).pivot(
            on="method", index="compound_id", values="prediction"
        )
        t = t.join(rt, on="compound_id", validate="1:1", maintain_order="left")
        assert t.height == 750
        qmaps = {"baseline_mean": p, "blend3_q55": q55}
        tmaps = {
            "baseline_mean": t["baseline:mean"].to_numpy(),
            "blend3_q55": t["blend3:q0.55"].to_numpy(),
        }
        neighbor = neighbor_oof.filter(
            (pl.col("endpoint") == ep) & (pl.col("method") == cfg["neighbor_method"])
        ).rename({"prediction": "neighbor"})
        nq = aligned_matrix(f, neighbor, ["neighbor"])[:, 0]
        nt = neighbor_test.filter(
            (pl.col("endpoint") == ep) & (pl.col("method") == cfg["neighbor_method"])
        ).select("compound_id", pl.col("prediction").alias("neighbor"))
        nt = (
            t.select("compound_id")
            .join(nt, on="compound_id", validate="1:1", maintain_order="left")["neighbor"]
            .to_numpy()
        )
        for strength in [0.1, 0.25]:
            name = f"neighbor_interval_blend_q55_{strength}"
            qmaps[name] = strength * nq + (1 - strength) * q55
            tmaps[name] = strength * nt + (1 - strength) * tmaps["blend3_q55"]
        for i, method in enumerate(methods):
            qmaps[method] = matrix[:, i]
            tmaps[method] = t[method].to_numpy()
            for strength in [0.25, 0.5]:
                name = method + f"_blend_q55_{strength}"
                qmaps[name] = strength * matrix[:, i] + (1 - strength) * q55
                tmaps[name] = strength * t[method].to_numpy() + (1 - strength) * tmaps["blend3_q55"]
        bootdata[ep] = dict(
            y=y,
            lo=lo,
            hi=hi,
            p=p,
            q55=q55,
            qmaps=qmaps,
            groups=np.array([groupmap[c] for c in f["compound_id"]]),
        )
        for name, q in qmaps.items():
            assert np.isfinite(q).all() and np.isfinite(tmaps[name]).all()
            for fold in [-1, 0, 1, 2, 3, 4]:
                ix = np.ones(len(y), dtype=bool) if fold == -1 else folds == fold
                score = official.rae_soft_threshold_absolute_error(
                    y[ix], q[ix], y_true_lower=lo[ix], y_true_upper=hi[ix]
                )
                rows.append(
                    dict(
                        endpoint=ep,
                        method=name,
                        fold=fold,
                        ST_RAE=float(score),
                        MAE=float(np.abs(y[ix] - q[ix]).mean()),
                        Spearman=float(spearmanr(y[ix], q[ix]).statistic),
                    )
                )
            preds.append(
                f.select("compound_id", "endpoint", "fold", "split_id", "y_true").with_columns(
                    pl.lit(name).alias("method"), pl.Series("prediction", q)
                )
            )
            testpreds.append(
                t.select("compound_id").with_columns(
                    pl.lit(ep).alias("endpoint"),
                    pl.lit(name).alias("method"),
                    pl.Series("prediction", tmaps[name]),
                    pl.lit(-1).alias("fold"),
                    pl.lit("test:five-fold-point-mean").alias("split_id"),
                )
            )
    score = pl.DataFrame(rows)
    summary = (
        score.filter(pl.col("fold") == -1)
        .group_by("method")
        .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
        .sort("ST_RAE")
    )
    allgroups = sorted(set(g for d in bootdata.values() for g in d["groups"]))
    samples = np.random.default_rng(cfg["seed"]).multinomial(
        len(allgroups), np.ones(len(allgroups)) / len(allgroups), size=1000
    )
    bs = {}
    for ep, d in bootdata.items():
        ww = samples[:, [allgroups.index(g) for g in d["groups"]]]
        mean = ww @ d["y"] / ww.sum(axis=1)
        den = (
            ww
            * (
                np.maximum(d["lo"][None, :] - mean[:, None], 0)
                + np.maximum(mean[:, None] - d["hi"][None, :], 0)
            )
        ).sum(axis=1)
        for method, q in d["qmaps"].items():
            for anchor in ["p", "q55"]:
                bs[ep, method, anchor] = (
                    ww @ (soft(q, d["lo"], d["hi"]) - soft(d[anchor], d["lo"], d["hi"])) / den
                )
    boot = []
    eps = list(bootdata)
    for method in qmaps:
        for anchor in ["p", "q55"]:
            for ep in eps + ["macro"]:
                draws = (
                    np.mean([bs[e, method, anchor] for e in eps], axis=0)
                    if ep == "macro"
                    else bs[ep, method, anchor]
                )
                low, high = np.quantile(draws, [0.025, 0.975])
                boot.append(
                    dict(
                        endpoint=ep,
                        method=method,
                        reference=anchor,
                        delta_low=float(low),
                        delta_high=float(high),
                    )
                )
    score.write_csv(out / "interval-metrics.csv")
    summary.write_csv(out / "interval-summary.csv")
    pl.DataFrame(boot).write_csv(out / "interval-bootstrap.csv")
    pl.concat(preds).write_parquet(out / "interval-comparison-oof.parquet")
    pl.concat(testpreds).write_parquet(out / "interval-comparison-test.parquet")
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            audit_passed=True,
            results=summary.to_dicts(),
            scope="Ordinary fixed UMAP OOF; method/blend selection exploratory. "
            "Bootstrap excludes selection/refitting uncertainty. No hidden-test scores.",
        ),
    )
    print(summary)


if __name__ == "__main__":
    main()
