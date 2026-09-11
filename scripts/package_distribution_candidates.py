"""Audit keyed ensemble predictions, register metrics and prepare local candidate files."""

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import psycopg

from cyp_challenge.external_training import snapshot_input
from cyp_challenge.storage import digest, write_json
from cyp_challenge.tracks import ENDPOINTS, validate


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    prediction_file = snapshot_input(
        inp, cfg.get("prediction_filename", "distribution-comparison-test.parquet")
    )
    p = pl.read_parquet(prediction_file)
    aliases = (
        pl.read_parquet(snapshot_input(inp, "aliases.parquet"))
        .filter(pl.col("role") == "test")
        .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
        .unique()
    )
    reference = snapshot_input(inp, "cyp-challenge-TEST-BLINDED.csv")
    ids = (
        pl.read_csv(reference)
        .select("SMILES", "Molecule_Name")
        .join(aliases, on="Molecule_Name", validate="1:1", maintain_order="left")
    )
    assert ids.height == 750 and ids["compound_id"].n_unique() == 750
    stats, checks = [], []
    for method, name in cfg["candidates"].items():
        part = p.filter(pl.col("method") == method)
        assert part.select("compound_id", "endpoint").n_unique() == 3000
        wide = part.pivot(on="endpoint", index="compound_id", values="prediction")
        file = out / (name + ".parquet")
        result = ids.join(wide, on="compound_id", validate="1:1", maintain_order="left").select(
            "SMILES", "Molecule_Name", *ENDPOINTS["regression"]
        )
        result.write_parquet(file)
        checks.append(
            dict(method=method, file=file.name, validation=validate(file, "regression", reference))
        )
    methods = cfg.get(
        "plot_methods",
        ["baseline:mean", "baseline:q0.55", "blend3:mean", "blend3:median", "blend3:q0.55"],
    )
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
    for ep, ax in zip(ENDPOINTS["regression"], axes.flat, strict=True):
        for method in methods:
            x = p.filter((pl.col("endpoint") == ep) & (pl.col("method") == method))[
                "prediction"
            ].to_numpy()
            assert len(x) == 750 and np.isfinite(x).all()
            stats.append(
                dict(
                    endpoint=ep,
                    method=method,
                    mean=float(x.mean()),
                    sd=float(x.std()),
                    below4=int((x < 4).sum()),
                    above6=int((x > 6).sum()),
                )
            )
            ax.step(np.sort(x), np.arange(1, len(x) + 1) / len(x), where="post", label=method)
        ax.set(title=ep.split("_")[0], xlabel="Predicted pIC50", ylabel="Test cumulative fraction")
        ax.legend(fontsize=8)
    fig.savefig(out / "test-distributions.png", dpi=150)
    plt.close(fig)
    pl.DataFrame(stats).write_csv(out / "test-distributions.csv")
    current = Path(os.environ["CYP_RUN_DIR"]).name
    count = 0
    with psycopg.connect(**cfg["database"]) as conn:
        status, manifest = conn.execute(
            "SELECT status,manifest FROM runs WHERE id=%s", (cfg["prediction_run"],)
        ).fetchone()
        assert status == "completed"
        assert digest(prediction_file) in {v["sha256"] for v in manifest["files"]}
        for name, parent in cfg["metrics"].items():
            source = snapshot_input(inp, name + ".csv")
            status, manifest = conn.execute(
                "SELECT status,manifest FROM runs WHERE id=%s", (parent,)
            ).fetchone()
            assert status == "completed" and digest(source) in {
                v["sha256"] for v in manifest["files"]
            }
            for r in pl.read_csv(source).iter_rows(named=True):
                for key in ["ST_RAE", "MAE", "Spearman"]:
                    conn.execute(
                        "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
                        (
                            current,
                            r["endpoint"],
                            name + ":repeat" + str(r.get("repeat", -1)),
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
            candidates=checks,
            registered_metrics=count,
            test_statistics=stats,
            scope="Local review candidates only; no remote submission or hidden-test score",
        ),
    )
    print(pl.DataFrame(stats))


if __name__ == "__main__":
    main()
