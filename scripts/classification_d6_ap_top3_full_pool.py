"""Compose a frozen AP-mass90 top-three proxy from complete-pool OOF atoms."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score

inputs = Path(os.environ["CYP_INPUT_DIR"])
outputs = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
summary = json.loads((inputs / cfg["ap_summary"]).read_text())
weights = summary["weights"]
members = cfg["members"]
selected = [member["name"] for member in members]
if len(set(selected)) != len(selected) or not all(
    name in weights and weights[name] > 0 for name in selected
):
    raise ValueError("The frozen AP member selection does not match the summary")
if sorted(selected, key=lambda name: -weights[name]) != selected:
    raise ValueError("Members are not in decreasing submitted-weight order")
captured_weight = sum(weights[name] for name in selected)
if not cfg["expected_weight_min"] < captured_weight < cfg["expected_weight_max"]:
    raise ValueError("Unexpected AP weight coverage")
normalized = np.asarray([weights[name] / captured_weight for name in selected])
outer = pl.read_parquet(inputs / cfg["outer_details"]).filter(pl.col("method") == "ap::mass90")
if outer.height != 5 or sorted(outer["held_fold"].to_list()) != list(range(5)):
    raise ValueError("Incomplete outer-excluded AP-mass90 weights")
outer_weights = {}
for row in outer.iter_rows(named=True):
    fold_weights = json.loads(row["weights"])
    values = np.asarray([fold_weights.get(name, 0.0) for name in selected])
    if values.sum() <= 0:
        raise ValueError(f"No selected top-three atom in held fold {row['held_fold']}")
    outer_weights[row["held_fold"]] = values / values.sum()

sources = {}
for member in members:
    key = member["source"]
    if key not in sources:
        sources[key] = pl.read_parquet(inputs / cfg[key])

output_rows = []
for strategy in cfg.get("strategies", ["umap", "scaffold"]):
    parts = []
    for member in members:
        part = (
            sources[member["source"]]
            .filter(
                (pl.col("split_id") == strategy)
                & (pl.col("model") == member["model"])
                & (pl.col("endpoint") == cfg["endpoint"])
            )
            .sort("compound_id")
        )
        if part.height != 6146 or part["compound_id"].n_unique() != 6146:
            raise ValueError(f"Incomplete atomic OOF: {strategy} {member['name']}")
        parts.append(part)
    reference = parts[0]
    for part in parts[1:]:
        if not (
            part["compound_id"].equals(reference["compound_id"])
            and part["fold"].equals(reference["fold"])
            and part["y_true"].equals(reference["y_true"])
        ):
            raise ValueError(f"Atomic OOF key, fold or label disagreement: {strategy}")
    folds = reference["fold"].to_numpy()
    rank_matrix = np.empty((reference.height, len(parts)))
    for column, part in enumerate(parts):
        values = part["score"].to_numpy()
        if not np.isfinite(values).all():
            raise ValueError("Nonfinite atomic score")
        for fold in np.unique(folds):
            mask = folds == fold
            rank_matrix[mask, column] = rankdata(values[mask], method="average") / (
                mask.sum() + 1.0
            )
    score = rank_matrix @ normalized
    output_rows.append(
        reference.select(
            "compound_id", "fold", "group", "y_true", "endpoint", "split_id"
        ).with_columns(pl.lit(cfg["model"]).alias("model"), pl.Series("score", score))
    )
    if strategy == "umap" or cfg.get("nested_selection") == "caruana_ap":
        nested = np.empty(reference.height)
        for fold in np.unique(folds):
            mask = folds == fold
            if cfg.get("nested_selection") == "caruana_ap":
                train = (~mask) & reference["y_true"].is_not_null().to_numpy()
                truth = reference["y_true"].to_numpy()[train].astype(int)
                values = rank_matrix[train]
                counts = np.zeros(len(members), dtype=int)
                running = np.zeros(len(truth))
                best_ap, best_weights = -np.inf, None
                for step in range(1, cfg["nested_steps"] + 1):
                    ap = np.asarray(
                        [
                            average_precision_score(truth, (running + values[:, index]) / step)
                            for index in range(len(members))
                        ]
                    )
                    choice = int(np.argmax(ap))
                    counts[choice] += 1
                    running += values[:, choice]
                    if ap[choice] > best_ap + 1e-12:
                        best_ap = ap[choice]
                        best_weights = counts.astype(float) / step
                nested[mask] = rank_matrix[mask] @ best_weights
                outer_weights[int(fold)] = best_weights
            else:
                nested[mask] = rank_matrix[mask] @ outer_weights[int(fold)]
        output_rows.append(
            reference.select(
                "compound_id", "fold", "group", "y_true", "endpoint", "split_id"
            ).with_columns(pl.lit(cfg["nested_model"]).alias("model"), pl.Series("score", nested))
        )

pl.concat(output_rows).write_parquet(outputs / "all-pool-oof.parquet")
(outputs / "summary.json").write_text(
    json.dumps(
        {
            "model": cfg["model"],
            "source_fullfit_weight_coverage": captured_weight,
            "normalized_member_weights": dict(zip(selected, normalized.tolist(), strict=True)),
            "outer_excluded_normalized_member_weights": {
                str(fold): dict(zip(selected, values.tolist(), strict=True))
                for fold, values in outer_weights.items()
            },
            "scope": cfg["scope"],
        },
        indent=2,
    )
    + "\n"
)
