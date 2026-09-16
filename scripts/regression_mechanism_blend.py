"""Audit fixed blends of CYP mechanism readouts with retained regression ensembles."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(index: int) -> Path:
    return next(inp.glob(f"{index:03d}-*"))


mechanism, latest, submitted, memberships, assignments = [
    pl.read_parquet(path(index)) for index in range(5)
]
spec = importlib.util.spec_from_file_location("official_scorer", path(5))
if spec is None or spec.loader is None:
    raise RuntimeError("Could not load the official scorer")
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)

keys = ["compound_id", "endpoint", "fold"]
latest = latest.sort(keys)
submitted = submitted.sort(keys)
if latest.height != submitted.height or latest.height != 6384:
    raise ValueError("Reference OOF coverage mismatch")
for column in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
    if not latest[column].equals(submitted[column]):
        raise ValueError(f"Reference OOF mismatch in {column}")
for (endpoint, fold), membership in memberships.group_by("endpoint", "outer_fold"):
    membership = membership.join(
        assignments.select("compound_id", "fold"), on="compound_id", validate="1:1"
    )
    expected = (membership["fold"].to_numpy() == fold) & ~membership["is_anchor"].to_numpy()
    if not np.array_equal(membership["validation"].to_numpy(), expected):
        raise ValueError(f"Parent membership mismatch for {endpoint} fold {fold}")

predictions = {
    "latest": latest["prediction"].to_numpy(),
    "submitted": submitted["prediction"].to_numpy(),
}
for method in sorted(mechanism["method"].unique()):
    if method == "control":
        continue
    frame = mechanism.filter(pl.col("method") == method).sort(keys)
    if frame.height != 6384 or frame.select(keys).is_duplicated().any():
        raise ValueError(f"Mechanism OOF coverage mismatch for {method}")
    for column in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
        if not latest[column].equals(frame[column]):
            raise ValueError(f"Mechanism OOF mismatch for {method} in {column}")
    prediction = frame["prediction"].to_numpy()
    predictions[method] = prediction
    for base in ["latest", "submitted"]:
        for weight in cfg["weights"]:
            name = f"{base}__{method}__{int(100 * weight)}pct"
            predictions[name] = (1 - weight) * predictions[base] + weight * prediction

aligned = latest.select(keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]).with_columns(
    [pl.Series(name, values) for name, values in predictions.items()]
)
rows = []
for (endpoint,), endpoint_frame in aligned.group_by("endpoint"):
    for subset in ["all", "neighbors"]:
        frame = (
            endpoint_frame
            if subset == "all"
            else endpoint_frame.filter(pl.col("is_neighbor"))
        )
        y = frame["y_true"].to_numpy()
        lower = frame["conf_low"].to_numpy()
        upper = frame["conf_high"].to_numpy()
        for method in predictions:
            prediction = frame[method].to_numpy()
            rows.append(
                {
                    "endpoint": endpoint,
                    "subset": subset,
                    "method": method,
                    "n": frame.height,
                    "ST-RAE": float(
                        scorer.rae_soft_threshold_absolute_error(
                            y,
                            prediction,
                            y_true_lower=lower,
                            y_true_upper=upper,
                        )
                    ),
                    "MAE": float(np.abs(y - prediction).mean()),
                    "Spearman": float(spearmanr(y, prediction).statistic),
                }
            )

metrics = pl.DataFrame(rows).sort("subset", "endpoint", "ST-RAE", "method")
macro = (
    metrics.group_by("method", "subset")
    .agg(pl.col("ST-RAE", "MAE", "Spearman").mean())
    .sort("subset", "ST-RAE", "method")
)
metrics.write_csv(out / "metrics.csv")
macro.write_csv(out / "macro.csv")
aligned.write_parquet(out / "aligned-oof.parquet")
(out / "verification.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "matched_per_method": 6384,
            "exact_keys_labels_bounds": True,
            "parent_membership_verified": True,
            "weight_search": False,
        },
        indent=2,
    )
    + "\n"
)
print(macro, flush=True)
