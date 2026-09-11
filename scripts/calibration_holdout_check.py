"""Fit calibrators inside each base-model held-out fold, split by UMAP clusters.

Both calibration and evaluation labels were excluded from that base model.
This trades calibration sample size and deployment fidelity for label isolation.
"""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from calibration_screen import apply_map, fit_map, soft
from rdkit import Chem
from scipy.stats import spearmanr

from cyp_challenge.external_training import snapshot_input
from cyp_challenge.storage import write_json


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def read(name):
        return pl.read_parquet(snapshot_input(inp, name + ".parquet"))

    data = read("regression-oof-diagnostics").join(
        read("importance-weights").select("compound_id", "endpoint", "weight"),
        on=["compound_id", "endpoint"],
        validate="1:1",
    )
    groups = dict(
        read("assignments")
        .filter(pl.col("strategy") == "umap")
        .select("compound_id", "group")
        .iter_rows()
    )
    keys = {
        c: Chem.MolToInchiKey(Chem.MolFromSmiles(s)).split("-")[0]
        for c, s in read("compounds").select("compound_id", "standardized_smiles").iter_rows()
    }
    spec = importlib.util.spec_from_file_location(
        "official", snapshot_input(inp, "custom_scoring_functions.py")
    )
    official = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(official)
    rows, predictions, parameters, membership = [], [], [], []
    for ep in sorted(data["endpoint"].unique()):
        f = data.filter(pl.col("endpoint") == ep).sort("compound_id")
        ids = f["compound_id"].to_numpy()
        p, y, lo, hi, w = [
            f[c].to_numpy().astype(float)
            for c in ["prediction", "y_true", "conf_low", "conf_high", "weight"]
        ]
        folds = f["fold"].to_numpy()
        group = np.array([groups[c] for c in ids])
        conn = np.array([keys[c] for c in ids])
        for repeat in range(cfg["repeats"]):
            partitions = []
            for fold in range(5):
                available = np.unique(group[folds == fold])
                shuffled = np.random.default_rng(cfg["seed"] + repeat * 100 + fold).permutation(
                    available
                )
                half = set(shuffled[: len(shuffled) // 2])
                side = np.array([g in half for g in group])
                for s in [False, True]:
                    va = (folds == fold) & (side == s)
                    tr = (folds == fold) & (side != s) & ~np.isin(conn, conn[va])
                    assert tr.sum() > 20 and va.sum() > 20
                    assert set(conn[tr]).isdisjoint(conn[va])
                    partitions.append((fold, s, tr, va))
                    membership.append(
                        f.filter(pl.Series(va))
                        .select("compound_id", "endpoint", "fold")
                        .with_columns(
                            pl.lit(repeat).alias("repeat"),
                            pl.lit(int(s)).alias("evaluation_half"),
                            pl.lit(int(tr.sum())).alias("calibration_n"),
                        )
                    )
            for name in cfg["methods"]:
                kind = name.removesuffix("_iw").removesuffix("_half")
                iw = name.endswith("_iw")
                q = np.full(len(p), np.nan)
                for fold, s, tr, va in partitions:
                    model = fit_map(
                        kind, p[tr], y[tr], lo[tr], hi[tr], w[tr] if iw else np.ones(tr.sum())
                    )
                    corrected = apply_map(model, p[va])
                    q[va] = 0.5 * (corrected + p[va]) if name.endswith("_half") else corrected
                    parameters.append(
                        dict(
                            endpoint=ep,
                            method=name,
                            repeat=repeat,
                            fold=fold,
                            half=int(s),
                            n=int(tr.sum()),
                            model=model,
                        )
                    )
                assert np.isfinite(q).all()
                for fold in [-1, 0, 1, 2, 3, 4]:
                    ix = np.ones(len(p), dtype=bool) if fold == -1 else folds == fold
                    score = float(
                        official.rae_soft_threshold_absolute_error(
                            y[ix], q[ix], y_true_lower=lo[ix], y_true_upper=hi[ix]
                        )
                    )
                    den = soft(np.full(ix.sum(), y[ix].mean()), lo[ix], hi[ix]).sum()
                    assert np.isclose(score, soft(q[ix], lo[ix], hi[ix]).sum() / den)
                    rows.append(
                        dict(
                            endpoint=ep,
                            method=name,
                            repeat=repeat,
                            fold=fold,
                            ST_RAE=score,
                            MAE=float(np.abs(q[ix] - y[ix]).mean()),
                            Spearman=float(spearmanr(y[ix], q[ix]).statistic),
                        )
                    )
                predictions.append(
                    f.select("compound_id", "endpoint", "fold", "split_id", "y_true").with_columns(
                        pl.Series("prediction", q),
                        pl.lit(name).alias("method"),
                        pl.lit(repeat).alias("repeat"),
                    )
                )
    metrics = pl.DataFrame(rows)
    macro = (
        metrics.filter(pl.col("fold") == -1)
        .group_by("method", "repeat")
        .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
        .sort("method", "repeat")
    )
    summary = (
        macro.group_by("method")
        .agg(
            pl.col("ST_RAE").mean(),
            pl.col("ST_RAE").min().alias("min"),
            pl.col("ST_RAE").max().alias("max"),
            pl.col("MAE", "Spearman").mean(),
        )
        .sort("ST_RAE")
    )
    metrics.write_csv(out / "holdout-metrics.csv")
    macro.write_csv(out / "holdout-macro-repeats.csv")
    summary.write_csv(out / "holdout-summary.csv")
    pl.concat(predictions).write_parquet(out / "holdout-predictions.parquet")
    pl.concat(membership).write_parquet(out / "holdout-membership.parquet")
    write_json(out / "parameters.json", parameters)
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            audit_passed=True,
            results=summary.to_dicts(),
            limitations=[
                "Calibration/evaluation are disjoint UMAP-cluster halves "
                "within each base-model holdout fold.",
                "Much smaller calibration sets than deployment; same-fold chemistry "
                "may be more representative than deployment calibration.",
                "Repeated splits reuse labels: range is stability, "
                "not an independent confidence interval.",
                "Method selection remains exploratory; no hidden test performance measured.",
            ],
        ),
    )
    print(summary)


if __name__ == "__main__":
    main()
