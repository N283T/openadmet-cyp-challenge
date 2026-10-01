"""Matched OOF prediction/residual correlations for the adopted trunk placements."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    selections = pl.read_parquet(inp / "selections.parquet")
    existing = pl.read_parquet(inp / "incumbent-oof.parquet")
    keys = ["compound_id", "endpoint", "fold"]
    rows, chosen_rows = [], []
    for target, setup in cfg["targets"].items():
        bank = pl.read_parquet(inp / setup["bank_oof"])
        for column, value in setup["filter"].items():
            bank = bank.filter(pl.col(column) == value)
        bank = bank.filter(pl.col("endpoint") == setup["endpoint"]).sort(keys)
        table = bank.select(*keys, "y_true", "conf_low", "conf_high")
        incumbent = existing.filter(
            (pl.col("target") == target) & (pl.col("scenario") == "current")
        )
        table = table.join(
            incumbent.select(*keys, pl.col("prediction").alias("incumbent")),
            on=keys,
            validate="1:1",
        ).sort(keys)
        references = {"incumbent_ensemble": table["incumbent"].to_numpy()}
        if target == "cyp1a2":
            specs = json.loads((inp / setup["bank_summary"]).read_text())["config"]["components"]
            for spec in specs:
                atom = pl.read_parquet(inp / f"{target}-{spec['name']}.parquet").filter(
                    pl.col("endpoint") == setup["endpoint"]
                )
                joined = table.select(*keys).join(
                    atom.select(*keys, "prediction"),
                    on=keys,
                    validate="1:1",
                    maintain_order="left",
                )
                references[spec["name"]] = joined["prediction"].to_numpy()
        else:
            base = pl.read_parquet(inp / "d6-base.parquet").filter(
                (pl.col("endpoint") == setup["endpoint"]) & (pl.col("scenario") == "current_two")
            )
            base = (
                table.select(*keys)
                .join(
                    base.select(*keys, "prediction"),
                    on=keys,
                    validate="1:1",
                    maintain_order="left",
                )["prediction"]
                .to_numpy()
            )
            references["submitted_two_member_core"] = base
            references["moment_matched_robust_branch"] = 2 * references["incumbent_ensemble"] - base
        trunk = pl.read_parquet(inp / setup["predictions"]).filter(pl.col("role") == "train")
        candidates = {}
        for name in sorted(trunk["model"].unique()):
            joined = table.select(*keys).join(
                trunk.filter(pl.col("model") == name).select(*keys, "prediction"),
                on=keys,
                validate="1:1",
                maintain_order="left",
            )
            candidates[name] = joined["prediction"].to_numpy()
        folds = table["fold"].to_numpy()
        selected = np.empty(table.height)
        selection = "placement" if target == "cyp1a2" else "bank_reselection"
        choices = selections.filter(
            (pl.col("target") == target) & (pl.col("selection") == selection)
        )
        if choices.height != 5:
            raise ValueError("Missing fold selections")
        for row in choices.iter_rows(named=True):
            model = row["chosen"].split("::", 1)[1]
            selected[folds == row["held_fold"]] = candidates[model][folds == row["held_fold"]]
            chosen_rows.append({"target": target, "fold": row["held_fold"], "model": model})
        candidates["selected_raw_trunk"] = selected
        y, low, high = (table[c].to_numpy() for c in ["y_true", "conf_low", "conf_high"])
        for model, p in candidates.items():
            if not np.isfinite(p).all():
                raise ValueError("Incomplete OOF predictions")
            for reference, q in references.items():
                if not np.isfinite(q).all():
                    raise ValueError("Incomplete reference OOF")
                for held in [-1, 0, 1, 2, 3, 4]:
                    mask = np.ones(table.height, dtype=bool) if held == -1 else folds == held
                    a, b, truth = p[mask], q[mask], y[mask]
                    interval_a = a - np.clip(a, low[mask], high[mask])
                    interval_b = b - np.clip(b, low[mask], high[mask])
                    rows.append(
                        {
                            "target": target,
                            "model": model,
                            "reference": reference,
                            "fold": held,
                            "n": int(mask.sum()),
                            "prediction_pearson": float(np.corrcoef(a, b)[0, 1]),
                            "prediction_spearman": float(spearmanr(a, b).statistic),
                            "residual_pearson": float(np.corrcoef(truth - a, truth - b)[0, 1]),
                            "signed_interval_residual_pearson": float(
                                np.corrcoef(interval_a, interval_b)[0, 1]
                            ),
                        }
                    )
    pl.DataFrame(rows).write_parquet(out / "correlations.parquet")
    pl.DataFrame(chosen_rows).write_parquet(out / "selected-models.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "residual_definition": "observed pIC50 minus prediction",
                "interval_residual_definition": (
                    "prediction minus its projection onto confidence interval"
                ),
                "selected_raw_trunk": "Fold-specific selected readout before ensemble mixing",
                "scope": (
                    "Descriptive fitted OOF correlation, conditional on prior recipe selection"
                ),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
