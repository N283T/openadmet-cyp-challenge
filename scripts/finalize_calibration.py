"""Audit saved corrections, create reviewable test variants, and register screen metrics."""

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import psycopg
from scipy.stats import spearmanr

from cyp_challenge.external_training import snapshot_input
from cyp_challenge.storage import digest, write_json


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def read(name):
        return pl.read_csv(snapshot_input(inp, name + ".csv"))

    t = pl.read_parquet(snapshot_input(inp, "calibration-test.parquet"))
    params = json.loads(snapshot_input(inp, "parameters.json").read_text())
    ts, transforms, outputs = [], [], []
    for ep in sorted(t["endpoint"].unique()):
        base = t.filter((pl.col("endpoint") == ep) & (pl.col("method") == "identity")).sort(
            "compound_id"
        )
        full = t.filter((pl.col("endpoint") == ep) & (pl.col("method") == "affine_st")).sort(
            "compound_id"
        )
        assert base["compound_id"].to_list() == full["compound_id"].to_list()
        m = next(
            r["model"]
            for r in params
            if r["endpoint"] == ep and r["method"] == "affine_st" and r["fold"] == -1
        )
        assert np.allclose(full["prediction"], m["a"] * base["prediction"] + m["b"])
        for strength in [0.0, 0.5, 1.0]:
            a, b = 1 + strength * (m["a"] - 1), strength * m["b"]
            pred = a * base["prediction"].to_numpy() + b
            assert np.isfinite(pred).all() and a > 0
            assert np.isclose(spearmanr(base["prediction"], pred).statistic, 1)
            name = {0.0: "identity", 0.5: "affine_st_half", 1.0: "affine_st"}[strength]
            transforms.append(dict(endpoint=ep, method=name, a=a, b=b))
            ts.append(
                dict(
                    endpoint=ep,
                    method=name,
                    mean=float(pred.mean()),
                    sd=float(pred.std()),
                    below4=int((pred < 4).sum()),
                )
            )
            outputs.append(
                base.with_columns(pl.Series("prediction", pred), pl.lit(name).alias("method"))
            )
    result = pl.concat(outputs)
    assert result.select("compound_id", "endpoint", "method").n_unique() == 9000
    result.write_parquet(out / "review-test-variants.parquet")
    pl.DataFrame(ts).write_csv(out / "test-distributions.csv")
    pl.DataFrame(transforms).write_csv(out / "transforms.csv")
    # Nonuniform scenario coordinates must be plotted at their actual locations.
    scenarios = read("sensitivity")
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
    for ep, ax in zip(sorted(scenarios["endpoint"].unique()), axes.flat, strict=True):
        f = scenarios.filter(pl.col("endpoint") == ep)
        shifts = sorted(f["shift"].unique())
        scales = sorted(f["scale"].unique())
        z = np.array(
            [
                [
                    f.filter((pl.col("shift") == s) & (pl.col("scale") == a))["ST_RAE"][0]
                    for a in scales
                ]
                for s in shifts
            ]
        )
        im = ax.pcolormesh(scales, shifts, z, shading="nearest")
        ax.scatter([1], [0], c="red", marker="x")
        ax.set(
            title=ep.split("_")[0] + ": OOF ST-RAE",
            xlabel="Scale around prediction mean",
            ylabel="Shift in pIC50",
        )
        fig.colorbar(im, ax=ax)
    fig.savefig(out / "affine-sensitivity.png", dpi=150)
    plt.close(fig)
    current = Path(os.environ["CYP_RUN_DIR"]).name
    count = 0
    with psycopg.connect(**cfg["database"]) as conn:
        for name, parent in [
            ("metrics", cfg["screen_run"]),
            ("holdout-metrics", cfg["holdout_run"]),
        ]:
            status, manifest = conn.execute(
                "SELECT status,manifest FROM runs WHERE id=%s", (parent,)
            ).fetchone()
            assert status == "completed"
            assert digest(snapshot_input(inp, name + ".csv")) in {
                v["sha256"] for v in manifest["files"]
            }
            for r in read(name).iter_rows(named=True):
                split = (
                    "umap:secondary-calibration"
                    if name == "metrics"
                    else "umap:within-holdout:repeat" + str(r["repeat"])
                )
                for key in ["ST_RAE", "MAE", "Spearman"]:
                    conn.execute(
                        "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
                        (
                            current,
                            r["endpoint"],
                            split,
                            str(r["fold"]),
                            r["method"] + ":" + key,
                            r[key],
                        ),
                    )
                    count += 1
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            audit_passed=True,
            registered_metrics=count,
            decision="Retain original baseline. Keep interval-affine and half-strength versions "
            "for prospective validation; no upload.",
            note="Review variants use full-OOF fitted maps. Their hidden-test score is unknown.",
        ),
    )
    print(pl.DataFrame(ts))
    print(pl.DataFrame(transforms))


if __name__ == "__main__":
    main()
