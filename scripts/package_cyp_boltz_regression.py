"""Deploy accepted trunk ensemble families, preserving best-submission moments."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from evaluate_cyp_boltz_incumbent import error, selected_weights
from evaluate_cyp_boltz_multi import replacements
from scipy.stats import spearmanr

from cyp_challenge.tracks import validate


def match(p, reference):
    if not np.isfinite(p).all() or p.std() <= 0:
        raise ValueError("Invalid prediction vector")
    return (p - p.mean()) / p.std() * reference.std() + reference.mean()


def align(table, frame, keys, column="prediction"):
    joined = table.select(keys).join(
        frame.select(*keys, column), on=keys, validate="1:1", maintain_order="left"
    )
    values = joined[column].to_numpy()
    if not np.isfinite(values).all():
        raise ValueError("Incomplete keyed predictions")
    return values


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    best = pl.read_csv(inp / "best.csv").sort("Molecule_Name")
    last = pl.read_csv(inp / "last.csv").sort("Molecule_Name")
    reference = pl.read_csv(inp / "reference.csv").sort("Molecule_Name")
    aliases = (
        pl.read_parquet(inp / "aliases.parquet")
        .filter(pl.col("role") == "test")
        .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
        .unique()
    )
    test = reference.join(aliases, on="Molecule_Name", validate="1:1", maintain_order="left")
    if (
        test.height != 750
        or test["compound_id"].null_count()
        or test["compound_id"].n_unique() != 750
    ):
        raise ValueError("Invalid test identity map")
    if not reference.equals(best.select(reference.columns)):
        raise ValueError("Best submission identity mismatch")
    candidate, decisions, comparisons, long_frames, checks = best.clone(), [], [], [], []
    keys = ["compound_id", "endpoint", "fold"]
    for target, setup in cfg["targets"].items():
        endpoint = setup["endpoint"]
        bank = pl.read_parquet(inp / setup["bank_oof"])
        for col, value in setup["filter"].items():
            bank = bank.filter(pl.col(col) == value)
        table = bank.filter(pl.col("endpoint") == endpoint).sort(keys)
        if table.height != setup["expected_rows"] or table.select(keys).is_duplicated().any():
            raise ValueError("Wrong OOF cohort")
        y, low, high, folds = (
            table[c].to_numpy() for c in ["y_true", "conf_low", "conf_high", "fold"]
        )
        current = table["prediction"].to_numpy()
        trunk = pl.read_parquet(inp / setup["predictions"])
        names = sorted(n for n in trunk["model"].unique() if not n.startswith("molecular_only"))
        z = np.column_stack(
            [
                align(
                    table, trunk.filter((pl.col("role") == "train") & (pl.col("model") == n)), keys
                )
                for n in names
            ]
        )

        def trunk_test(name, bank=trunk):
            part = bank.filter((pl.col("role") == "test") & (pl.col("model") == name))
            if part.height != 3750 or sorted(part["fold"].unique()) != list(range(5)):
                raise ValueError("Missing five-fold test predictions")
            return np.column_stack(
                [align(test, part.filter(pl.col("fold") == f), ["compound_id"]) for f in range(5)]
            ).mean(axis=1)

        old = best[endpoint].to_numpy()
        if target == "cyp1a2":
            specs = json.loads((inp / setup["bank_summary"]).read_text())["config"]["components"]
            atoms = [s["name"] for s in specs]
            x = np.column_stack(
                [
                    align(
                        table,
                        pl.read_parquet(inp / f"{target}-{name}.parquet").filter(
                            pl.col("endpoint") == endpoint
                        ),
                        keys,
                    )
                    for name in atoms
                ]
            )
            w = np.full(x.shape, 1 / len(atoms))
            if np.max(np.abs(x.mean(axis=1) - current)) > 1e-9:
                raise ValueError("OOF reconstruction mismatch")
            blind_atoms = pl.read_parquet(inp / "atomic-test.parquet").filter(
                pl.col("endpoint") == endpoint
            )
            xt = np.column_stack(
                [
                    align(test, blind_atoms.filter(pl.col("member") == n), ["compound_id"])
                    for n in atoms
                ]
            )
            baseline_raw = xt.mean(axis=1)
            scores = []
            for cap in [0, 1, 2, 3]:
                p = np.empty(len(y))
                for held in range(5):
                    q, _ = replacements(z, x, w, current, low, high, folds != held, cap)
                    p[folds == held] = q[folds == held]
                scores.append((float(error(low, high, p).mean()), cap))
            _, cap = min(scores)
            _, choices = replacements(z, x, w, current, low, high, np.ones(len(y), dtype=bool), cap)
            raw = baseline_raw.copy()
            for j, k in choices:
                raw += (trunk_test(names[k]) - xt[:, j]) / len(atoms)
            decisions.append(
                {
                    "target": target,
                    "method": "multiple_replaced",
                    "cap": cap,
                    "inner_scores": scores,
                    "replacements": [
                        {"old": atoms[j], "new": names[k], "weight": 0.1} for j, k in choices
                    ],
                }
            )
        else:
            base = pl.read_parquet(inp / "d6-base.parquet").filter(
                (pl.col("endpoint") == endpoint) & (pl.col("scenario") == "current_two")
            )
            core = align(table, base, keys)
            x = np.column_stack([core, 2 * current - core])
            k = min(
                range(len(names)), key=lambda k: (float(error(low, high, z[:, k]).mean()), names[k])
            )
            extension = np.column_stack([x, z[:, k]])
            w = selected_weights(extension, low, high, folds, "caruana")
            base_test = pl.read_parquet(inp / "d6-core-test.parquet").filter(
                pl.col("endpoint") == endpoint
            )
            core_test = align(test, base_test, ["compound_id"])
            robust_test = align(
                test,
                pl.read_parquet(inp / "d6-robust-test.parquet"),
                ["compound_id"],
                "raw_prediction",
            )
            robust_test = match(robust_test, core_test)
            baseline_raw = 0.5 * core_test + 0.5 * robust_test
            raw = np.column_stack([core_test, robust_test, trunk_test(names[k])]) @ w
            decisions.append(
                {
                    "target": target,
                    "method": "single_trunk_bank_extension",
                    "trunk": names[k],
                    "weights": w.tolist(),
                    "members": [
                        "submitted_two_member_core",
                        "moment_matched_robust_branch",
                        names[k],
                    ],
                }
            )
        baseline = match(baseline_raw, old)
        difference = float(np.max(np.abs(baseline - old)))
        if difference > 1e-8:
            raise ValueError(f"Blind incumbent reconstruction failed: {target} {difference}")
        checks.append({"target": target, "baseline_max_absolute_difference": difference})
        placed = match(raw, old)
        candidate = candidate.with_columns(pl.Series(endpoint, placed))
        long_frames.append(
            test.select("compound_id", "Molecule_Name").with_columns(
                pl.lit(endpoint).alias("endpoint"),
                pl.lit("blind_five_fold_mean").alias("split_id"),
                pl.lit(-1).alias("fold"),
                pl.Series("raw_prediction", raw),
                pl.Series("prediction", placed),
            )
        )
    for endpoint in best.columns[2:]:
        p, b, previous = (f[endpoint].to_numpy() for f in [candidate, best, last])
        comparisons.append(
            {
                "endpoint": endpoint,
                "changed_vs_best": int(np.count_nonzero(p != b)),
                "changed_vs_last_sent": int(np.count_nonzero(p != previous)),
                "mean": float(p.mean()),
                "sample_sd": float(p.std(ddof=1)),
                "mean_absolute_delta_vs_last": float(np.abs(p - previous).mean()),
                "spearman_vs_best": float(spearmanr(p, b).statistic),
            }
        )
        if endpoint not in [v["endpoint"] for v in cfg["targets"].values()] and not candidate[
            endpoint
        ].equals(best[endpoint]):
            raise ValueError("Changed rollback target")
    candidate.write_csv(out / "predictions.csv")
    candidate.write_parquet(out / "predictions.parquet")
    pl.concat(long_frames).write_parquet(out / "changed-target-test-predictions.parquet")
    pl.DataFrame(comparisons).write_csv(out / "target-comparison.csv")
    validation = validate(out / "predictions.csv", "regression", inp / "reference.csv")
    if not validation["valid"]:
        raise ValueError(validation)
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "decisions": decisions,
                "baseline_checks": checks,
                "comparison": comparisons,
                "validation": validation,
                "deployment": (
                    "Select composition on all development OOF; average five trunk test folds. "
                    "Original atoms use retained full-fit predictions. Match final moments "
                    "to the best submitted endpoint, without new calibration search. "
                    "Development selection is conditional, not independent validation."
                ),
            },
            indent=2,
        )
        + "\n"
    )
    print(
        json.dumps({"decisions": decisions, "checks": checks, "comparison": comparisons}),
        flush=True,
    )


if __name__ == "__main__":
    main()
