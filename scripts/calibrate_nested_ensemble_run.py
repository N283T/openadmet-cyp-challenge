"""Inner-only affine calibration of saved nested Caruana ensembles; no model fits."""

import importlib.util
import json
import os
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl
import psycopg
from scipy.optimize import linprog
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
keys = ["compound_id", "endpoint", "fold"]
outer = pl.read_parquet(inp / cfg["outer_oof"]).sort(keys)
jeremy = pl.read_parquet(inp / cfg["jeremy_oof"]).sort(keys)
assert outer.height == jeremy.height == 6384
for name in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
    assert outer[name].equals(jeremy[name]), name
assert outer.select(keys).unique().height == outer.height
rows, parameters, reconstructed = [], [], []


def score(frame, prediction):
    return float(
        scorer.rae_soft_threshold_absolute_error(
            frame["y_true"].to_numpy(),
            prediction,
            y_true_lower=frame["conf_low"].to_numpy(),
            y_true_upper=frame["conf_high"].to_numpy(),
        )
    )


def fit(inner, prediction, variant):
    # ST-RAE denominator depends only on inner y. Minimizing summed distances to
    # confidence intervals is therefore exactly minimizing the official objective.
    n = len(prediction)
    identity = np.eye(n)
    upper = np.column_stack([prediction, np.ones(n), -identity])
    lower = np.column_stack([-prediction, -np.ones(n), -identity])
    constraints = np.concatenate([upper, lower])
    bounds_rhs = np.concatenate([inner["conf_high"].to_numpy(), -inner["conf_low"].to_numpy()])
    scale_bounds = (1.0, 1.0) if variant == "offset" else tuple(cfg["scale_bounds"])
    bounds = [scale_bounds, tuple(cfg["offset_bounds"])] + [(0.0, None)] * n
    result = linprog(
        np.r_[0.0, 0.0, np.ones(n)],
        A_ub=constraints,
        b_ub=bounds_rhs,
        bounds=bounds,
        method="highs",
        options={"dual_feasibility_tolerance": 1e-9, "primal_feasibility_tolerance": 1e-9},
    )
    assert result.success, result.message
    scale, offset = map(float, result.x[:2])
    prediction_after = scale * prediction + offset
    before, after = score(inner, prediction), score(inner, prediction_after)
    assert after <= before + 1e-8
    summed_error = np.clip(prediction_after - inner["conf_high"].to_numpy(), 0, None).sum()
    summed_error += np.clip(inner["conf_low"].to_numpy() - prediction_after, 0, None).sum()
    assert np.isclose(summed_error, result.fun, atol=1e-7)
    return dict(
        scale=scale,
        offset=offset,
        inner_before=before,
        inner_after=after,
        success=bool(result.success),
        status=int(result.status),
        message=result.message,
        iterations=int(result.nit),
        summed_interval_error=float(result.fun),
    ), prediction_after


for item in cfg["groups"]:
    inner = pl.read_parquet(inp / item["inner"])
    member = pl.read_parquet(inp / item["membership"])
    selected = json.loads((inp / item["selected"]).read_text())
    history = json.loads((inp / item["history"]).read_text())
    ep, fold = selected["endpoint"], selected["fold"]
    assert member["endpoint"].unique().to_list() == [ep]
    assert member["outer_fold"].unique().to_list() == [fold]
    assert member["compound_id"].n_unique() == member.height
    holdout = member.filter(pl.col("role") == "outer_validation")
    train = member.filter(pl.col("role") == "inner_train")
    validation = member.filter(pl.col("role") == "inner_validation")
    assert set(inner["compound_id"]) == set(validation["compound_id"])
    assert not set(inner["compound_id"]) & set(holdout["compound_id"])
    assert not set(train["group"]) & set(validation["group"])
    assert not validation["is_anchor"].any()
    assert (
        member.filter((pl.col("fold") == fold) & ~pl.col("is_anchor"))
        .select("compound_id")
        .sort("compound_id")
        .equals(holdout.select("compound_id").sort("compound_id"))
    )
    target = outer.filter((pl.col("endpoint") == ep) & (pl.col("fold") == fold))
    assert set(target["compound_id"]) == set(holdout["compound_id"])
    counts = Counter(step["selected"] for step in history[: selected["best_step"]])
    assert counts == {candidate["name"]: candidate["count"] for candidate in selected["selected"]}
    assert sum(counts.values()) == selected["best_step"]
    prediction = np.zeros(inner.height)
    for candidate in selected["selected"]:
        assert np.isclose(candidate["weight"], candidate["count"] / selected["best_step"])
        prediction += candidate["weight"] * inner[candidate["name"]].to_numpy()
    assert np.isclose(score(inner, prediction), selected["inner_score"], rtol=0, atol=1e-12)
    assert np.isclose(
        history[selected["best_step"] - 1]["score"], selected["inner_score"], rtol=0, atol=1e-12
    )
    saved = inner.select(keys + ["y_true", "conf_low", "conf_high"]).with_columns(
        pl.lit(fold).alias("outer_fold"), pl.Series("reconstructed_prediction", prediction)
    )
    for variant in ["offset", "affine"]:
        param, calibrated_inner = fit(inner, prediction, variant)
        parameters.append(
            dict(
                endpoint=ep,
                fold=fold,
                variant=variant,
                inner_rows=inner.height,
                outer_rows=target.height,
                selected_step=selected["best_step"],
                **param,
            )
        )
        calibrated_outer = param["scale"] * target["prediction"].to_numpy() + param["offset"]
        assert np.isfinite(calibrated_outer).all()
        rows.append(
            target.with_columns(
                pl.Series("prediction", calibrated_outer),
                pl.lit("nested_caruana24_" + variant).alias("method"),
            )
        )
        saved = saved.with_columns(pl.Series(variant + "_prediction", calibrated_inner))
    reconstructed.append(saved)
pl.DataFrame(parameters).write_csv(out / "calibration-parameters.csv")
(out / "calibration-parameters.json").write_text(json.dumps(parameters, indent=2))
pl.concat(reconstructed).write_parquet(out / "reconstructed-inner-predictions.parquet")
variants = pl.concat(rows)
for variant in ["offset", "affine"]:
    prediction = variants.filter(pl.col("method") == "nested_caruana24_" + variant).sort(keys)
    assert prediction.height == 6384 and prediction.select(keys).equals(outer.select(keys))
    prediction.write_parquet(out / (variant + "-oof.parquet"))
base = outer.with_columns(
    outer["prediction"].alias("uncalibrated"), jeremy["prediction"].alias("jeremy_frozen_tabpfn")
)
for variant in ["offset", "affine"]:
    prediction = variants.filter(pl.col("method") == "nested_caruana24_" + variant).sort(keys)
    base = base.with_columns(prediction["prediction"].alias(variant))
for variant in ["uncalibrated", "offset", "affine"]:
    base = base.with_columns(
        pl.mean_horizontal(variant, "jeremy_frozen_tabpfn").alias(variant + "_jeremy_equal")
    )
methods = [
    "uncalibrated",
    "offset",
    "affine",
    "jeremy_frozen_tabpfn",
    "uncalibrated_jeremy_equal",
    "offset_jeremy_equal",
    "affine_jeremy_equal",
]
metrics = []
for (ep,), group in base.group_by("endpoint"):
    for subset in ["all", "neighbors"]:
        selected = group if subset == "all" else group.filter(pl.col("is_neighbor"))
        for fold in [-1, 0, 1, 2, 3, 4]:
            h = selected if fold == -1 else selected.filter(pl.col("fold") == fold)
            for method in methods:
                p = h[method].to_numpy()
                metrics.append(
                    dict(
                        endpoint=ep,
                        subset=subset,
                        fold=fold,
                        n=h.height,
                        method=method,
                        ST_RAE=score(h, p),
                        MAE=float(np.abs(h["y_true"].to_numpy() - p).mean()),
                        Spearman=float(spearmanr(h["y_true"].to_numpy(), p).statistic),
                    )
                )
metrics = pl.DataFrame(metrics)
metrics.write_csv(out / "metrics.csv")
macro = (
    metrics.filter(pl.col("fold") == -1)
    .group_by("method", "subset")
    .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
    .sort("subset", "ST_RAE")
)
macro.write_csv(out / "macro.csv")
base.write_parquet(out / "aligned-oof.parquet")
with psycopg.connect(**cfg["database"]) as conn:
    for row in metrics.iter_rows(named=True):
        for metric in ["ST_RAE", "MAE", "Spearman"]:
            conn.execute(
                "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
                (
                    Path(os.environ["CYP_RUN_DIR"]).name,
                    row["endpoint"],
                    "nested_calibration:" + row["subset"],
                    str(row["fold"]),
                    row["method"] + ":" + metric,
                    row[metric],
                ),
            )
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            calibration_fits=len(parameters),
            variants=["offset", "affine"],
            rows_per_variant=6384,
            original_selections_unchanged=True,
            calibration_fit_data="Saved inner validation ensemble only; no outer targets",
            score_reconstruction_tolerance=1e-12,
            optimization="Exact bounded convex interval-loss linear program",
            outer_evaluation_did_not_choose_variant=True,
        ),
        indent=2,
    )
)
print(macro, flush=True)
