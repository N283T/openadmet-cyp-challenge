"""Independent integrity, fold, and bootstrap audit of selected D6 screen rows."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
keys = ["compound_id", "endpoint", "fold"]


def read(index: int) -> pl.DataFrame:
    paths = list(inp.glob(f"{index:03d}-*"))
    if len(paths) != 1:
        raise ValueError(f"Expected one input at index {index}, found {paths}")
    return pl.read_parquet(paths[0])


def st_rae(frame: pl.DataFrame, prediction: np.ndarray) -> float:
    y = frame["y_true"].to_numpy()
    lower = frame["conf_low"].to_numpy()
    upper = frame["conf_high"].to_numpy()
    error = np.maximum(lower - prediction, 0) + np.maximum(prediction - upper, 0)
    center = np.full_like(y, y.mean())
    denominator = (np.maximum(lower - center, 0) + np.maximum(center - upper, 0)).sum()
    if denominator <= 0:
        raise ValueError("Non-positive ST-RAE denominator")
    return float(error.sum() / denominator)


oof, latest, current = [read(index) for index in range(3)]
latest = latest.filter(pl.col("endpoint") == cfg["endpoint"]).sort(keys)
current = current.filter(
    (pl.col("endpoint") == cfg["endpoint"]) & (pl.col("component_count") == 2)
).sort(keys)
if latest.height != current.height or latest.height != cfg["expected_rows"]:
    raise ValueError("Baseline coverage mismatch")
for column in keys + ["y_true", "conf_low", "conf_high"]:
    if not latest[column].equals(current[column]):
        raise ValueError(f"Baseline mismatch in {column}")

expected_methods = {tuple(method) for method in cfg["methods"]}
observed_methods = set(zip(oof["view"], oof["method"], strict=True))
if observed_methods != expected_methods:
    raise ValueError(f"Method inventory mismatch: {observed_methods ^ expected_methods}")
if oof.height != cfg["expected_rows"] * len(expected_methods):
    raise ValueError("Screen OOF row count mismatch")
for (view, method), frame in oof.group_by("view", "method"):
    frame = frame.sort(keys)
    for column in keys + ["y_true", "conf_low", "conf_high"]:
        if not frame[column].equals(latest[column]):
            raise ValueError(f"OOF mismatch in {view}/{method}/{column}")

base_prediction = latest["prediction"].to_numpy()
y = latest["y_true"].to_numpy()
baseline_score = st_rae(latest, base_prediction)
baseline_spearman = float(spearmanr(y, base_prediction).statistic)
rng = np.random.default_rng(cfg["bootstrap_seed"])
audit_rows = []
fold_rows = []
bootstrap_rows = []
candidate_frames = []
for candidate in cfg["candidates"]:
    frame = oof.filter(
        (pl.col("view") == candidate["view"]) & (pl.col("method") == candidate["method"])
    ).sort(keys)
    prediction = frame["prediction"].to_numpy()
    if candidate["placement"] == "match_latest_moments":
        prediction = (prediction - prediction.mean()) / prediction.std() * base_prediction.std()
        prediction += base_prediction.mean()
    elif candidate["placement"] != "raw":
        raise ValueError(candidate["placement"])
    blended = (1 - candidate["weight"]) * base_prediction + candidate["weight"] * prediction
    name = candidate["name"]
    candidate_frames.append(
        frame.select(*keys, "group", "is_neighbor", "y_true", "conf_low", "conf_high").with_columns(
            pl.lit(name).alias("candidate"), pl.Series("prediction", blended)
        )
    )
    candidate_score = st_rae(latest, blended)
    candidate_spearman = float(spearmanr(y, blended).statistic)
    audit_rows.append(
        {
            **candidate,
            "ST_RAE": candidate_score,
            "delta_ST_RAE": candidate_score - baseline_score,
            "Spearman": candidate_spearman,
            "delta_Spearman": candidate_spearman - baseline_spearman,
            "prediction_correlation": float(np.corrcoef(prediction, base_prediction)[0, 1]),
            "residual_correlation": float(np.corrcoef(y - prediction, y - base_prediction)[0, 1]),
        }
    )
    for fold in range(5):
        mask = latest["fold"].to_numpy() == fold
        fold_frame = latest.filter(pl.Series(mask))
        fold_rows.append(
            {
                "candidate": name,
                "fold": fold,
                "n": int(mask.sum()),
                "baseline_ST_RAE": st_rae(fold_frame, base_prediction[mask]),
                "candidate_ST_RAE": st_rae(fold_frame, blended[mask]),
            }
        )

    group_to_indices = {
        group: indices["index"].to_numpy()
        for (group,), indices in frame.with_row_index("index").group_by("group")
    }
    groups = np.asarray(sorted(group_to_indices))
    deltas = []
    for _ in range(cfg["bootstrap_replicates"]):
        sampled_groups = rng.choice(groups, size=len(groups), replace=True)
        indices = np.concatenate([group_to_indices[group] for group in sampled_groups])
        sampled = latest[indices]
        deltas.append(st_rae(sampled, blended[indices]) - st_rae(sampled, base_prediction[indices]))
    values = np.asarray(deltas)
    bootstrap_rows.append(
        {
            "candidate": name,
            "replicates": cfg["bootstrap_replicates"],
            "delta_mean": float(values.mean()),
            "delta_q025": float(np.quantile(values, 0.025)),
            "delta_q50": float(np.quantile(values, 0.5)),
            "delta_q975": float(np.quantile(values, 0.975)),
            "probability_improvement": float(np.mean(values < 0)),
        }
    )

audit = pl.DataFrame(audit_rows).sort("ST_RAE")
folds = (
    pl.DataFrame(fold_rows)
    .with_columns((pl.col("candidate_ST_RAE") - pl.col("baseline_ST_RAE")).alias("delta_ST_RAE"))
    .sort("candidate", "fold")
)
bootstrap = pl.DataFrame(bootstrap_rows).sort("delta_q50")
audit.write_csv(out / "audit.csv")
folds.write_csv(out / "folds.csv")
bootstrap.write_csv(out / "bootstrap.csv")
pl.concat(candidate_frames).write_parquet(out / "candidate-oof.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "baseline_ST_RAE": baseline_score,
            "baseline_Spearman": baseline_spearman,
            "current_two_component_ST_RAE": st_rae(current, current["prediction"].to_numpy()),
            "integrity": "Exact keys, labels, bounds, method inventory, and row counts verified.",
            "inference_scope": (
                "Bootstrap is conditional on saved OOF predictions and post-screen candidate "
                "selection; intervals are not selection-adjusted."
            ),
        },
        indent=2,
    )
    + "\n"
)
print(audit, flush=True)
print(folds, flush=True)
print(bootstrap, flush=True)
