"""Target-wise exhaustive pruning of the submitted five-component regression mean."""

import importlib.util
import itertools
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
keys = ["compound_id", "endpoint", "fold"]
truth = ["y_true", "conf_low", "conf_high", "is_neighbor"]

spec = importlib.util.spec_from_file_location("official_scorer", inp / cfg["scorer"])
assert spec is not None and spec.loader is not None
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)

frames = []
for component in cfg["components"]:
    frame = pl.read_parquet(inp / component["file"]).sort(keys)
    assert frame.height == 6384
    assert not frame.select(keys).is_duplicated().any()
    frames.append(frame)

base = frames[0]
for frame in frames[1:]:
    for column in keys + truth:
        assert frame[column].equals(base[column]), column

names = [component["name"] for component in cfg["components"]]
predictions = np.column_stack(
    [
        frame[component["column"]].to_numpy()
        for frame, component in zip(frames, cfg["components"], strict=True)
    ]
)
assert np.isfinite(predictions).all()


def score(frame: pl.DataFrame, prediction: np.ndarray) -> float:
    return float(
        scorer.rae_soft_threshold_absolute_error(
            frame["y_true"].to_numpy(),
            prediction,
            y_true_lower=frame["conf_low"].to_numpy(),
            y_true_upper=frame["conf_high"].to_numpy(),
        )
    )


subset_rows = []
all_subset_rows = []
selection_rows = []
crossfit_rows = []
crossfit_predictions = []
rng = np.random.default_rng(cfg["bootstrap_seed"])

indexed = base.with_row_index("row_index")
for (endpoint,), group in indexed.group_by("endpoint", maintain_order=True):
    group = group.sort(keys)
    ix = group["row_index"].to_numpy()
    x = predictions[ix]
    folds = group["fold"].to_numpy()
    combinations = {
        size: list(itertools.combinations(range(len(names)), size))
        for size in range(1, len(names) + 1)
    }

    for size, subsets in combinations.items():
        full_scores = [(score(group, x[:, subset].mean(axis=1)), subset) for subset in subsets]
        all_subset_rows.extend(
            {
                "endpoint": endpoint,
                "component_count": size,
                "components": "+".join(names[i] for i in subset),
                "ST_RAE": value,
            }
            for value, subset in full_scores
        )
        full_score, full_subset = min(full_scores, key=lambda item: (item[0], item[1]))
        subset_rows.append(
            {
                "endpoint": endpoint,
                "component_count": size,
                "selection_scope": "full_oof_fit",
                "components": "+".join(names[i] for i in full_subset),
                "ST_RAE": full_score,
            }
        )

        p_crossfit = np.empty(group.height, dtype=float)
        for fold in sorted(set(folds)):
            train = folds != fold
            valid = ~train
            train_frame = group.filter(pl.Series(train))
            choices = [
                (score(train_frame, x[train][:, subset].mean(axis=1)), subset) for subset in subsets
            ]
            train_score, selected = min(choices, key=lambda item: (item[0], item[1]))
            p_crossfit[valid] = x[valid][:, selected].mean(axis=1)
            valid_frame = group.filter(pl.Series(valid))
            selection_rows.append(
                {
                    "endpoint": endpoint,
                    "component_count": size,
                    "outer_fold": int(fold),
                    "components": "+".join(names[i] for i in selected),
                    "selection_ST_RAE": train_score,
                    "outer_ST_RAE": score(valid_frame, p_crossfit[valid]),
                }
            )

        cf_score = score(group, p_crossfit)
        crossfit_rows.append(
            {
                "endpoint": endpoint,
                "component_count": size,
                "ST_RAE": cf_score,
            }
        )
        crossfit_predictions.append(
            group.select(keys + truth).with_columns(
                pl.Series("prediction", p_crossfit),
                pl.lit(size).alias("component_count"),
            )
        )

    endpoint_cf = {
        row["component_count"]: row["ST_RAE"]
        for row in crossfit_rows
        if row["endpoint"] == endpoint
    }
    best_size = min(endpoint_cf, key=lambda size: (endpoint_cf[size], size))
    pred_by_size = {
        frame["component_count"][0]: frame["prediction"].to_numpy()
        for frame in crossfit_predictions
        if frame["endpoint"][0] == endpoint
    }
    best_pred = pred_by_size[best_size]
    for size in range(1, len(names) + 1):
        deltas = []
        for _ in range(cfg["bootstrap_replicates"]):
            sample = rng.integers(0, group.height, group.height)
            sampled = group[sample]
            deltas.append(
                score(sampled, pred_by_size[size][sample]) - score(sampled, best_pred[sample])
            )
        values = np.asarray(deltas)
        for row in crossfit_rows:
            if row["endpoint"] == endpoint and row["component_count"] == size:
                row.update(
                    {
                        "best_component_count": best_size,
                        "delta_vs_best": endpoint_cf[size] - endpoint_cf[best_size],
                        "delta_ci_low": float(np.quantile(values, 0.025)),
                        "delta_ci_high": float(np.quantile(values, 0.975)),
                    }
                )

all_subset_frame = pl.DataFrame(all_subset_rows)
shared_full = (
    all_subset_frame.group_by("component_count", "components")
    .agg(pl.col("ST_RAE").mean().alias("macro_ST_RAE"))
    .sort("component_count", "macro_ST_RAE", "components")
)
shared_best_rows = []
shared_fold_rows = []
all_folds = base["fold"].to_numpy()
for size in range(1, len(names) + 1):
    candidates = list(itertools.combinations(range(len(names)), size))
    full_best = shared_full.filter(pl.col("component_count") == size).row(0, named=True)
    p_crossfit = np.empty(base.height, dtype=float)
    for fold in sorted(set(all_folds)):
        train = all_folds != fold
        valid = ~train
        train_frame = base.filter(pl.Series(train)).with_row_index("local_index")
        endpoint_indices = {
            endpoint: group["local_index"].to_numpy()
            for (endpoint,), group in train_frame.group_by("endpoint", maintain_order=True)
        }
        scored = []
        for subset in candidates:
            component_prediction = predictions[train][:, subset].mean(axis=1)
            endpoint_scores = []
            for _endpoint, local_indices in endpoint_indices.items():
                endpoint_frame = train_frame[local_indices]
                endpoint_scores.append(score(endpoint_frame, component_prediction[local_indices]))
            scored.append((float(np.mean(endpoint_scores)), subset))
        train_score, selected = min(scored, key=lambda item: (item[0], item[1]))
        p_crossfit[valid] = predictions[valid][:, selected].mean(axis=1)
        shared_fold_rows.append(
            {
                "component_count": size,
                "outer_fold": int(fold),
                "components": "+".join(names[i] for i in selected),
                "selection_macro_ST_RAE": train_score,
            }
        )
    endpoint_scores = []
    indexed_base = base.with_row_index("row_index")
    for (_endpoint,), group in indexed_base.group_by("endpoint", maintain_order=True):
        row_indices = group["row_index"].to_numpy()
        endpoint_scores.append(score(group, p_crossfit[row_indices]))
    shared_best_rows.append(
        {
            "component_count": size,
            "crossfit_macro_ST_RAE": float(np.mean(endpoint_scores)),
            "full_oof_components": full_best["components"],
            "full_oof_macro_ST_RAE": full_best["macro_ST_RAE"],
        }
    )

pl.DataFrame(subset_rows).sort("endpoint", "component_count").write_csv(out / "full-oof-best.csv")
pl.DataFrame(all_subset_rows).sort("endpoint", "component_count", "ST_RAE").write_csv(
    out / "all-subsets.csv"
)
pl.DataFrame(selection_rows).sort("endpoint", "component_count", "outer_fold").write_csv(
    out / "fold-selections.csv"
)
pl.DataFrame(crossfit_rows).sort("endpoint", "component_count").write_csv(
    out / "crossfit-scores.csv"
)
pl.DataFrame(shared_best_rows).sort("component_count").write_csv(out / "shared-scores.csv")
pl.DataFrame(shared_fold_rows).sort("component_count", "outer_fold").write_csv(
    out / "shared-fold-selections.csv"
)
crossfit_frame = pl.concat(crossfit_predictions).sort(
    "endpoint", "component_count", "fold", "compound_id"
)
crossfit_frame.write_parquet(out / "crossfit-predictions.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "component_names": names,
            "selection": (
                "For each endpoint and size, select the best equal-weight subset on the "
                "other four folds and evaluate the held-out fold."
            ),
            "full_oof_note": (
                "Full-OOF best subsets are deployment-fit diagnostics, not held-out estimates."
            ),
            "bootstrap": {
                "replicates": cfg["bootstrap_replicates"],
                "seed": cfg["bootstrap_seed"],
                "unit": "endpoint-labelled compound row",
                "scope": "conditional on saved predictions and the cross-fitted selection results",
            },
        },
        indent=2,
    )
)

print(pl.DataFrame(crossfit_rows).sort("endpoint", "component_count"), flush=True)
print(pl.DataFrame(subset_rows).sort("endpoint", "component_count"), flush=True)
