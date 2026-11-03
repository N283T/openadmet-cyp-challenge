"""Super-ensemble of per-checkpoint, per-depth CheMeleon readouts for one regression endpoint.

Every recorded depth screen (`launch_openadmet_depth_screen.py`) scored one OpenADMET CheMeleon
checkpoint at one message-passing depth with a LightGBM and a TabPFN readout on the incumbent's own
OOF rows and folds. This script collects all of them and compares, with the incumbent as reference:
the members alone, simple and per-axis means, nested greedy selection (free and anchored on the
incumbent), and fixed-weight blends into the incumbent. Scores are pooled ST-RAE, and also ST-RAE on
the analog-proxy rows that stand in for the analog-expansion test set. Nothing is uploaded.
"""

import json
import os
import re
from pathlib import Path

import numpy as np
import polars as pl
import regression_bank_selection as sel

KEYS = ["compound_id", "fold"]
READOUTS = {"emb--all--lgbm_l1": "lgbm", "emb--variance:360--tabpfn": "tabpfn"}
PURPOSE = re.compile(r"^Screen (\S+) CheMeleon depth (\d+) readouts")


def collect(root, endpoint):
    """Latest screen run per (checkpoint, depth) with its OOF and blind-mean tables."""
    latest = {}
    for run in sorted(root.iterdir()):
        intent, config = run / "intent.json", run / "config.json"
        if not (intent.exists() and config.exists() and (run / "result.json").exists()):
            continue
        match = PURPOSE.match(json.loads(intent.read_text()).get("purpose", ""))
        if not match or json.loads((run / "result.json").read_text())["status"] != "completed":
            continue
        if not json.loads(config.read_text())["method"].startswith("openadmet_depth_member_screen"):
            continue
        latest[(match.group(1), int(match.group(2)))] = run
    members = {}
    for (name, depth), run in sorted(latest.items()):
        oof = pl.read_parquet(run / "outputs/predictions/oof.parquet").filter(
            pl.col("endpoint") == endpoint
        )
        blind = pl.read_parquet(run / "outputs/predictions/blind-mean.parquet").filter(
            pl.col("endpoint") == endpoint
        )
        for model, readout in READOUTS.items():
            members[f"{name}:{depth}:{readout}"] = (
                oof.filter(pl.col("model") == model),
                blind.filter(pl.col("model") == model),
                run.name,
            )
    return members


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    endpoint = cfg["endpoint"]
    control = cfg["control"]
    frame = pl.read_parquet(inp / control["file"])
    for column, value in control["filter"].items():
        frame = frame.filter(pl.col(column) == value)
    frame = frame.filter(pl.col("endpoint") == endpoint).sort("compound_id")
    assert frame.height == control["rows"]
    table = frame.select(
        *KEYS, "y_true", "conf_low", "conf_high", pl.col("prediction").alias("inc")
    )
    members = collect(Path(cfg["runs_root"]), endpoint)
    names = list(members)
    for name in names:
        part = members[name][0].select(*KEYS, pl.col("prediction").alias(name))
        table = table.join(part, on=KEYS, how="left", validate="1:1")
    table = table.sort("compound_id")
    assert table.null_count().sum_horizontal().item() == 0, "member rows do not cover the cohort"
    y, low, high, folds = (table[c].to_numpy() for c in ("y_true", "conf_low", "conf_high", "fold"))
    incumbent = table["inc"].to_numpy()
    matrix = table.select(names).to_numpy()
    analog = set(pl.read_parquet(inp / "analog-pairs.parquet")["analog_id"].to_list())
    is_analog = np.array([i in analog for i in table["compound_id"].to_list()])
    print(len(names), "members;", int(is_analog.sum()), "analog rows", flush=True)

    def scores(prediction):
        return dict(
            st_rae=sel.st_rae(y, low, high, prediction),
            analog=sel.st_rae(y[is_analog], low[is_analog], high[is_analog], prediction[is_analog]),
            other=sel.st_rae(
                y[~is_analog], low[~is_analog], high[~is_analog], prediction[~is_analog]
            ),
            folds_better=int(
                sum(
                    sel.st_rae(
                        y[folds == f], low[folds == f], high[folds == f], prediction[folds == f]
                    )
                    < sel.st_rae(
                        y[folds == f], low[folds == f], high[folds == f], incumbent[folds == f]
                    )
                    for f in sorted(set(folds))
                )
            ),
        )

    rows = [dict(procedure="incumbent", **scores(incumbent))]
    member_rows = []
    for j, name in enumerate(names):
        member_rows.append(dict(member=name, run=members[name][2], **scores(matrix[:, j])))
    pl.DataFrame(member_rows).write_csv(out / "members.csv")

    def column_mask(predicate):
        return np.array([predicate(n.split(":")) for n in names])

    means = {"mean_all": np.ones(len(names), bool)}
    for readout in ("lgbm", "tabpfn"):
        means[f"mean_{readout}"] = column_mask(lambda p, r=readout: p[2] == r)
    for depth in range(1, 7):
        means[f"mean_depth{depth}"] = column_mask(lambda p, d=depth: int(p[1]) == d)
    for checkpoint in sorted({n.split(":")[0] for n in names}):
        means[f"mean_ckpt_{checkpoint}"] = column_mask(lambda p, c=checkpoint: p[0] == c)
    for label, mask in means.items():
        rows.append(dict(procedure=label, **scores(matrix[:, mask].mean(axis=1))))

    full = np.column_stack([incumbent, matrix])
    nested = {}
    for procedure in ("greedy_free", "greedy_anchored"):
        blend = np.empty(len(y))
        for outer in sorted(set(folds)):
            fit, held = folds != outer, folds == outer
            weights, steps = sel.procedure_weights(
                full[fit], y[fit], low[fit], high[fit], folds[fit], procedure, cfg
            )
            blend[held] = full[held] @ weights
        nested[procedure] = blend
        rows.append(dict(procedure=f"nested_{procedure}", **scores(blend)))
    # Members only: selection that never sees the incumbent, i.e. how far the pool reaches alone.
    alone = {"members_best_single": np.empty(len(y)), "members_greedy": np.empty(len(y))}
    alone.update({f"members_top{k}_mean": np.empty(len(y)) for k in cfg["k_grid"] if k > 1})
    chosen = []
    for outer in sorted(set(folds)):
        fit, held = folds != outer, folds == outer
        fit_scores = [
            sel.st_rae(y[fit], low[fit], high[fit], matrix[fit, j]) for j in range(len(names))
        ]
        order = np.argsort(fit_scores)
        alone["members_best_single"][held] = matrix[held, order[0]]
        weights, steps = sel.procedure_weights(
            matrix[fit], y[fit], low[fit], high[fit], folds[fit], "greedy_free", cfg
        )
        alone["members_greedy"][held] = matrix[held] @ weights
        chosen.append(
            dict(
                outer=int(outer),
                steps=int(steps),
                best_single=names[order[0]],
                members={names[j]: round(float(w), 3) for j, w in enumerate(weights) if w > 0},
            )
        )
        for k in cfg["k_grid"]:
            if k > 1:
                alone[f"members_top{k}_mean"][held] = matrix[held][:, order[:k]].mean(axis=1)
    for label, blend in alone.items():
        rows.append(dict(procedure=f"nested_{label}", **scores(blend)))
    (out / "members-only-selection.json").write_text(json.dumps(chosen, indent=2) + "\n")
    weights, steps = sel.procedure_weights(full, y, low, high, folds, "greedy_anchored", cfg)
    final = {
        ("incumbent" if j == 0 else names[j - 1]): round(float(w), 4)
        for j, w in enumerate(weights)
        if w > 0
    }
    mean_all = matrix.mean(axis=1)
    for weight in cfg["blend_weights"]:
        rows.append(
            dict(
                procedure=f"incumbent+{weight}*mean_all",
                **scores((1 - weight) * incumbent + weight * mean_all),
            )
        )
    pl.DataFrame(rows).write_csv(out / "procedures.csv")
    (out / "final-anchored-weights.json").write_text(json.dumps(final, indent=2) + "\n")
    print(pl.DataFrame(rows).sort("st_rae"))
    blind = None
    for name in names:
        part = members[name][1].select("compound_id", pl.col("prediction").alias(name))
        blind = part if blind is None else blind.join(part, on="compound_id", how="inner")
    blind.write_parquet(out / "member-blind-predictions.parquet")


if __name__ == "__main__":
    main()
