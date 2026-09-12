"""Verify inner-only HPO decisions and compare prespecified three-model means."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
keys = ["compound_id", "endpoint", "fold"]
base = None
for name, filename in cfg["models"].items():
    frame = pl.read_parquet(inp / filename).sort(keys)
    assert frame.height == frame.select(keys).n_unique() == 6384
    columns = keys + ["y_true", "conf_low", "conf_high", "is_neighbor", "split_id"]
    if base is None:
        base = frame.select(columns)
    assert base.select(columns).equals(frame.select(columns))
    base = base.with_columns(frame["prediction"].alias(name))
base = base.with_columns(
    pl.mean_horizontal("caruana", "frozen_pfn", "fixed_lgbm").alias("fixed_three"),
    pl.mean_horizontal("caruana", "frozen_pfn", "tuned_lgbm").alias("tuned_three"),
)
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
rows = []
for (endpoint,), group in base.group_by("endpoint"):
    for subset in ["all", "neighbors"]:
        h = group if subset == "all" else group.filter(pl.col("is_neighbor"))
        for fold in [-1, 0, 1, 2, 3, 4]:
            g = h if fold == -1 else h.filter(pl.col("fold") == fold)
            for method in list(cfg["models"]) + ["fixed_three", "tuned_three"]:
                rows.append(
                    dict(
                        endpoint=endpoint,
                        subset=subset,
                        fold=fold,
                        method=method,
                        ST_RAE=float(
                            scorer.rae_soft_threshold_absolute_error(
                                g["y_true"].to_numpy(),
                                g[method].to_numpy(),
                                y_true_lower=g["conf_low"].to_numpy(),
                                y_true_upper=g["conf_high"].to_numpy(),
                            )
                        ),
                    )
                )
metrics = pl.DataFrame(rows)
metrics.write_csv(out / "metrics.csv")
macro = (
    metrics.filter(pl.col("fold") == -1)
    .group_by("method", "subset")
    .agg(pl.col("ST_RAE").mean())
    .sort("subset", "ST_RAE")
)
macro.write_csv(out / "macro.csv")
base.write_parquet(out / "aligned-oof.parquet")
base.select(keys + ["y_true", "conf_low", "conf_high", "is_neighbor", "split_id"]).with_columns(
    base["tuned_three"].alias("prediction"),
    pl.lit("caruana_jeremy_pfn_tuned_ftlgbm_equal").alias("method"),
).write_parquet(out / "external-oof.parquet")
membership = pl.read_parquet(inp / cfg["membership"])
reference = pl.read_parquet(inp / cfg["reference_membership"])
assert membership.sort(keys + ["outer_fold"]).equals(reference.sort(keys + ["outer_fold"]))
summary = json.loads((inp / cfg["summary"]).read_text())
encoder_checks = json.loads((inp / cfg["encoder_exclusion_audit"]).read_text())
assert len(encoder_checks) == 5 and all(c["outer_labels_excluded"] for c in encoder_checks)
assert sum(c["outer_compounds"] for c in encoder_checks) == 6076
selections = {(r["fold"], r["endpoint"]): r for r in summary["selections"]}
audit = []
for record in cfg["trial_files"]:
    trials = [json.loads(line) for line in (inp / record["file"]).read_text().splitlines()]
    assert len(trials) == 24 and [t["trial"] for t in trials] == list(range(24))
    assert np.isfinite([t["score"] for t in trials]).all()
    winner = min(trials, key=lambda t: t["score"])
    actual = selections[(record["fold"], record["endpoint"])]
    assert actual["trial"] == winner["trial"] and actual["params"] == winner["params"]
    assert actual["score"] == winner["score"]
    audit.append(
        dict(
            endpoint=record["endpoint"],
            fold=record["fold"],
            trial=winner["trial"],
            objective=winner["params"]["objective"],
            alpha=winner["params"].get("alpha"),
            baseline_inner_score=trials[0]["score"],
            winner_inner_score=winner["score"],
        )
    )
pl.DataFrame(audit, infer_schema_length=None).write_csv(out / "inner-selections.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            matched_oof=6384,
            inner_memberships_identical_to_caruana=True,
            finite_trials=480,
            selections_are_inner_argmin=True,
            weights_fitted_on_oof=False,
            inner_encoder_label_reuse=True,
        ),
        indent=2,
    )
)
lines = [
    "# Jeremy FT embedding + log2fc Rustuna readout",
    "",
    f"Training run `{cfg['fit_run']}` evaluates 24 LightGBM trials per outer fold and endpoint. "
    "The grouped 20% inner holdout matches the existing Caruana membership exactly; anchor groups "
    "remain in inner training. The saved winner is the minimum official inner ST-RAE in all "
    "20 fits, and all 480 trial scores are finite. Winners refit on all outer training labels.",
    "",
    "The FT encoder previously used labels from the readout's inner holdout, so inner HPO "
    "scores can be optimistic. This is not nested encoder selection. Each encoder's outer "
    "holdout was verified against parent assignments; those outer labels were excluded from "
    "encoder fitting, checkpoint selection and readout HPO. Upstream pretraining still lacks "
    "a molecule manifest, so this remains an upstream-model diagnostic.",
    "",
    "The comparison retains the fixed LGBM readout and uses the same prescribed equal "
    "three-way mean: original 24-candidate Caruana + frozen Jeremy TabPFN + FT/log2fc LGBM. "
    "No ensemble weights were fitted to OOF predictions.",
    "",
    "| Method | All ST-RAE | Neighbors ST-RAE |",
    "|---|---:|---:|",
]
for method in ["fixed_lgbm", "tuned_lgbm", "fixed_three", "tuned_three"]:
    values = {r["subset"]: r["ST_RAE"] for r in macro.filter(pl.col("method") == method).to_dicts()}
    lines.append(f"| {method} | {values['all']:.6f} | {values['neighbors']:.6f} |")
lines += [
    "",
    "| Endpoint | Fixed LGBM | Tuned LGBM | Fixed three | Tuned three |",
    "|---|---:|---:|---:|---:|",
]
for ep in sorted(metrics["endpoint"].unique()):
    v = {
        r["method"]: r["ST_RAE"]
        for r in metrics.filter(
            (pl.col("endpoint") == ep) & (pl.col("fold") == -1) & (pl.col("subset") == "all")
        ).to_dicts()
    }
    lines.append(
        f"| {ep.split('_')[0]} | {v['fixed_lgbm']:.6f} | {v['tuned_lgbm']:.6f} | "
        f"{v['fixed_three']:.6f} | {v['tuned_three']:.6f} |"
    )
(out / "report.md").write_text("\n".join(lines) + "\n")
print(macro, flush=True)
