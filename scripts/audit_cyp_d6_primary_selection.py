"""Backtest primary-conditioned label selection without inventing missing outcomes."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from cyp_primary_supervision import load_primary, propensity
from review_cyp_d6_causal_probe import calls_at_rate


def projected(y, calls, prior, weights):
    tpr = np.average(calls[y == 1], weights=weights[y == 1])
    fpr = np.average(calls[y == 0], weights=weights[y == 0])
    tp, fn, fp, tn = prior * tpr, prior * (1 - tpr), (1 - prior) * fpr, (1 - prior) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    compounds = pl.read_parquet(inp / "compounds.parquet").sort("compound_id")
    aliases = pl.read_parquet(inp / "aliases.parquet")
    measured = pl.read_parquet(inp / "measurements.parquet")
    table = compounds.join(
        measured.filter(pl.col("endpoint") == "CYP2D6_is_TDI").select(
            "compound_id", pl.col("value").alias("label")
        ),
        on="compound_id",
        how="left",
    )
    primary, _ = load_primary(inp / "primary-screen.csv", aliases, table)
    blind_ids = aliases.filter(pl.col("role") == "test")["compound_id"].unique()
    development = np.flatnonzero(~table["compound_id"].is_in(blind_ids).to_numpy())
    blind = np.flatnonzero(table["compound_id"].is_in(blind_ids).to_numpy())
    assert not np.isfinite(primary[blind]).any()
    y = table["label"].to_numpy()
    primary_d6 = primary[:, 0]
    bins = np.where(np.isfinite(primary_d6), np.digitize(primary_d6, [-2, -1, -0.5, 0]), 5)
    source = table.select("compound_id").with_columns(
        pl.Series("label", y),
        pl.Series("primary_log2fc", primary_d6),
        pl.Series("primary_bin", bins),
        pl.Series("is_development", np.isin(np.arange(len(table)), development)),
    )
    source.filter(pl.col("is_development")).write_parquet(out / "selection-source.parquet")
    source.filter(pl.col("is_development")).group_by("primary_bin").agg(
        pl.len().alias("pool_n"),
        pl.col("label").is_finite().sum().alias("labelled_n"),
        (pl.col("label") == 1).sum().alias("positive_n"),
    ).sort("primary_bin").write_csv(out / "selection-coverage-by-primary.csv")
    assignments = pl.read_parquet(inp / "assignments.parquet")
    weight_tables = {}
    for split in ["umap", "scaffold"]:
        ordered = table.select("compound_id").join(
            assignments.filter(pl.col("strategy") == split).select("compound_id", "fold", "group"),
            on="compound_id",
            how="left",
            maintain_order="left",
        )
        folds = ordered["fold"].to_numpy()
        prop = np.full(len(table), np.nan)
        for fold in range(5):
            fit = development[folds[development] != fold]
            valid = development[folds[development] == fold]
            assert set(ordered["group"][fit]).isdisjoint(ordered["group"][valid])
            prop[valid] = propensity(primary_d6[fit], np.isfinite(y[fit]), primary_d6[valid])
        weights = 1 / np.maximum(prop, 1 / cfg["weight_cap"])
        weight_table = ordered.with_columns(
            pl.Series("label", y),
            pl.Series("propensity", prop),
            pl.Series("ipw", weights),
            pl.Series("primary_bin", bins),
        )
        weight_tables[split] = weight_table
        weight_table.filter(
            pl.col("compound_id").is_in(table[development]["compound_id"])
        ).write_parquet(out / f"selection-propensity-{split}.parquet")
    comparisons, subgroup_rows, swaps = [], [], []
    for case in cfg["cases"]:
        parts = []
        for item in [case["before"], case["after"]]:
            raw = pl.read_parquet(inp / item["file"])
            model_col = "method" if "method" in raw.columns else "model"
            score_col = "prediction" if "prediction" in raw.columns else "score"
            part = raw.filter(pl.col(model_col) == item["model"]).select(
                "compound_id", "fold", "group", "y_true", pl.col(score_col).alias("score")
            )
            part = part.join(
                weight_tables[case["split"]].select("compound_id", "ipw", "primary_bin"),
                on="compound_id",
                how="left",
                validate="1:1",
            ).sort("compound_id")
            assert part["ipw"].is_finite().all()
            parts.append(part)
        a, b = parts
        assert a["compound_id"].equals(b["compound_id"]) and a["fold"].equals(b["fold"])
        np.testing.assert_allclose(a["y_true"], b["y_true"], equal_nan=True)
        ca = calls_at_rate(a["score"].to_numpy(), a["fold"].to_numpy(), case["before"]["rate"])
        cb = calls_at_rate(b["score"].to_numpy(), b["fold"].to_numpy(), case["after"]["rate"])
        observed = np.isfinite(a["y_true"].to_numpy())
        truth = a["y_true"].to_numpy()[observed].astype(int)
        groups = a["group"].to_numpy()[observed]
        for weighting in ["uniform", "selection_ipw"]:
            w = np.ones(len(truth)) if weighting == "uniform" else a["ipw"].to_numpy()[observed]
            indices_by_group = [np.flatnonzero(groups == g) for g in np.unique(groups)]
            rng = np.random.default_rng(20261002)
            deltas = []
            for _ in range(cfg["bootstrap_samples"]):
                draw = rng.integers(0, len(indices_by_group), len(indices_by_group))
                ix = np.concatenate([indices_by_group[i] for i in draw])
                if len(np.unique(truth[ix])) != 2:
                    continue
                deltas.append(
                    projected(truth[ix], cb[observed][ix], 0.07, w[ix])
                    - projected(truth[ix], ca[observed][ix], 0.07, w[ix])
                )
            before = projected(truth, ca[observed], 0.07, w)
            after = projected(truth, cb[observed], 0.07, w)
            comparisons.append(
                dict(
                    case=case["name"],
                    split_id=case["split"],
                    scope=case["scope"],
                    weighting=weighting,
                    lb_delta=case["lb_delta"],
                    oof_delta=after - before,
                    concordant=bool(np.sign(after - before) == np.sign(case["lb_delta"])),
                    ci_low=float(np.quantile(deltas, 0.025)),
                    ci_high=float(np.quantile(deltas, 0.975)),
                    effective_labelled_n=float(w.sum() ** 2 / (w * w).sum()),
                    assumptions="MAR given measured primary response; hypothetical prior.07; "
                    "conditional cached fits",
                )
            )
        for value in range(6):
            subset = observed & (a["primary_bin"].to_numpy() == value)
            truth_part = a["y_true"].to_numpy()[subset]
            adds = cb[subset] & ~ca[subset]
            removals = ca[subset] & ~cb[subset]
            subgroup_rows.append(
                dict(
                    case=case["name"],
                    primary_bin=value,
                    labelled=int(subset.sum()),
                    positives=int(truth_part.sum()),
                    added_calls=int(adds.sum()),
                    removed_calls=int(removals.sum()),
                    added_positive=int(((truth_part == 1) & adds).sum()),
                    removed_positive=int(((truth_part == 1) & removals).sum()),
                )
            )
        for i in np.flatnonzero(ca != cb):
            swaps.append(
                dict(
                    case=case["name"],
                    compound_id=a["compound_id"][int(i)],
                    fold=int(a["fold"][int(i)]),
                    split_id=case["split"],
                    endpoint="CYP2D6_is_TDI",
                    label=a["y_true"][int(i)],
                    primary_bin=int(a["primary_bin"][int(i)]),
                    before_call=bool(ca[i]),
                    after_call=bool(cb[i]),
                )
            )
    pl.DataFrame(comparisons).write_csv(out / "historical-concordance.csv")
    pl.DataFrame(subgroup_rows).write_csv(out / "primary-conditioned-swap-yields.csv")
    pl.DataFrame(swaps).write_parquet(out / "keyed-historical-swaps.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "comparisons": comparisons,
                "limits": "No weighting recovers unobserved outcomes; label selection can "
                "depend on unmeasured outcomes, violating MAR. Historical successes have "
                "partial-pool OOF and failures complete-pool OOF, explicitly separate scopes. "
                "No blind individual labels are inferred, no candidate is created.",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
