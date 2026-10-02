"""Select the diagnosed D6 AP-block procedure for a blind-call review only."""

import json
import os
import runpy
from pathlib import Path

import numpy as np
import polars as pl
from evaluate_cyp_boltz_tdi_diagnostics import compose, greedy_block, quality, rank_fold


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    scratch = Path(os.environ["TMPDIR"]) / "reconstruction"
    scratch.mkdir()
    os.environ["CYP_OUTPUT_DIR"] = str(scratch)
    ns = runpy.run_path(
        str(Path(__file__).with_name("evaluate_classification_depth_whole_members.py"))
    )
    os.environ["CYP_OUTPUT_DIR"] = str(out)
    y, folds = ns["d6_y"], ns["d6_folds"]
    ids = ns["d6_saved"].select("compound_id", "fold")
    current = ns["d6_current"]
    raw = pl.read_parquet(inp / cfg["trunk_predictions"]["d6"]).filter(pl.col("role") == "train")
    names = sorted(n for n in raw["model"].unique() if not n.startswith("molecular_only"))
    z = np.column_stack(
        [
            rank_fold(
                ids.join(
                    raw.filter(pl.col("model") == name).select("compound_id", "fold", "prediction"),
                    on=["compound_id", "fold"],
                    validate="1:1",
                    maintain_order="left",
                )["prediction"].to_numpy(),
                folds,
            )
            for name in names
        ]
    )
    if not np.isfinite(z).all():
        raise ValueError("Incomplete labelled predictions")
    rate, prior = cfg["d6"]["call_rate"], cfg["d6"]["target_prior"]
    base_rank = rank_fold(current, folds)
    options = {("current", 0, 0.0): current.copy()}
    for count in cfg["block_sizes"]:
        for kind, weight in [
            ("block_add", 0.0),
            ("individual_add", 0.0),
            *[("blend", w) for w in cfg["blend_weights"]],
        ]:
            options[(kind, count, weight)] = np.empty(len(y))
    fold_choices = []
    for held in sorted(np.unique(folds)):
        train, valid = folds != held, folds == held
        for count in cfg["block_sizes"]:
            block, chosen = greedy_block(z, names, y, folds, train, count, "AP", rate, prior)
            values = {
                ("block_add", count, 0.0): compose(ns, "d6", [block]),
                ("individual_add", count, 0.0): compose(ns, "d6", [z[:, i] for i in chosen]),
            }
            for w in cfg["blend_weights"]:
                values[("blend", count, w)] = (1 - w) * base_rank + w * block
            for key, score in values.items():
                options[key][valid] = score[valid]
            fold_choices.append(
                {
                    "held_fold": int(held),
                    "count": count,
                    "models": json.dumps([names[i] for i in chosen]),
                }
            )
    mask = np.ones(len(y), dtype=bool)
    metrics = [
        {
            "kind": key[0],
            "count": key[1],
            "weight": key[2],
            "inner_AP": quality(y, score, folds, mask, "AP", rate, prior),
        }
        for key, score in options.items()
    ]
    best = max(
        options,
        key=lambda key: (
            quality(y, options[key], folds, mask, "AP", rate, prior),
            key[0] == "current",
            -key[1],
            -key[2],
        ),
    )
    kind, count, weight = best
    chosen_names = []
    if count:
        _, indices = greedy_block(z, names, y, folds, mask, count, "AP", rate, prior)
        chosen_names = [names[i] for i in indices]
    decision = {
        "kind": kind,
        "count": count,
        "weight": weight,
        "models": chosen_names,
        "inner_AP": quality(y, options[best], folds, mask, "AP", rate, prior),
        "diagnostic_source": cfg["diagnostic_source"],
        "positive_calls": cfg["positive_calls"],
        "rate_policy": (
            "Keep historical 49/750 to isolate score-recipe changes; no rate optimization"
        ),
        "test_aggregation": (
            "Average five raw fold predictions per readout, rank across 750, "
            "then average selected readout ranks"
        ),
        "scope": (
            "Full-development selection of the previously diagnosed AP-block procedure; "
            "selection score is descriptive, not independent validation; no upload authorized"
        ),
    }
    pl.DataFrame(metrics).write_parquet(out / "inner-selection-metrics.parquet")
    pl.DataFrame(fold_choices).write_parquet(out / "inner-block-choices.parquet")
    (out / "decision.json").write_text(json.dumps(decision, indent=2) + "\n")
    print(json.dumps(decision, indent=2), flush=True)


if __name__ == "__main__":
    main()
