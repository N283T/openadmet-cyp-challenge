"""Compress the CYP3A4 Caruana component and evaluate it on outer-fold OOF data."""

import importlib.util
import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoint = cfg["endpoint"]


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


assignments = read("assignments").sort("compound_id")
labels = read("measurements")
compounds = read("compounds")
ids = assignments.select("compound_id")


def matrix(name: str) -> np.ndarray:
    table = read(name)
    columns = sorted(c for c in table.columns if c != "compound_id")
    x = (
        ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
        .select(columns)
        .to_numpy()
        .astype(np.float32)
    )
    assert np.isfinite(x).all()
    return x


xc, xm, xt = matrix("class_emb"), matrix("mixed_emb"), matrix("tree")
fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
smiles = ids.join(compounds, on="compound_id", validate="1:1", maintain_order="left")[
    "standardized_smiles"
]
xf = np.stack([fpgen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in smiles])
index = {compound_id: i for i, compound_id in enumerate(assignments["compound_id"])}

target = (
    labels.filter(pl.col("endpoint") == endpoint)
    .join(assignments, on="compound_id", validate="1:1")
    .sort("compound_id")
)
target_indices = np.array([index[c] for c in target["compound_id"]])

spec = importlib.util.spec_from_file_location("official_scorer", inp / cfg["scorer"])
assert spec is not None and spec.loader is not None
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)


def score(frame: pl.DataFrame, prediction: np.ndarray) -> float:
    return float(
        scorer.rae_soft_threshold_absolute_error(
            frame["y_true"].to_numpy(),
            prediction,
            y_true_lower=frame["conf_low"].to_numpy(),
            y_true_upper=frame["conf_high"].to_numpy(),
        )
    )


component_frames = []
for item in cfg["components"]:
    frame = pl.read_parquet(inp / item["file"]).filter(pl.col("endpoint") == endpoint)
    component_frames.append(frame.sort("fold", "compound_id"))
base = component_frames[0]
keys = ["compound_id", "endpoint", "fold"]
truth = ["y_true", "conf_low", "conf_high", "is_neighbor"]
for frame in component_frames[1:]:
    for column in keys + truth:
        assert frame[column].equals(base[column]), column
component_by_name = {
    item["name"]: frame for item, frame in zip(cfg["components"], component_frames, strict=True)
}
ensemble_names = cfg.get("ensemble_components", list(component_by_name))
assert ensemble_names[0] == cfg["components"][0]["name"]
assert set(ensemble_names).issubset(component_by_name)
ensemble_frames = [component_by_name[name] for name in ensemble_names]
other_sum = sum(frame["prediction"].to_numpy() for frame in ensemble_frames[1:])
original_caruana = base["prediction"].to_numpy()
original_ensemble = (original_caruana + other_sum) / len(ensemble_frames)

fold_records = {int(item["fold"]): item for item in cfg["folds"]}
assert set(fold_records) == set(range(5))
selection_by_fold = {}
candidate_specs = {}
saved_predictions = {}
reconstruction_rows = []
fixed_atomic_rows = {"class_quantile_0.65": [], "morgan_quantile_0.65": []}

for fold in range(5):
    record = fold_records[fold]
    selection = json.loads((inp / record["selection"]).read_text())
    selected = sorted(selection["selected"], key=lambda item: (-item["weight"], item["name"]))
    selection_by_fold[fold] = selected
    fold_rows = base.filter(pl.col("fold") == fold)
    fold_target = target.filter((pl.col("fold") == fold) & ~pl.col("is_anchor"))
    assert fold_rows["compound_id"].equals(fold_target["compound_id"])
    ix = np.array([index[c] for c in fold_rows["compound_id"]])
    nested = pl.read_parquet(inp / record["nested"])
    xn = (
        ids.join(
            nested.filter(pl.col("outer_fold") == fold),
            on="compound_id",
            validate="1:1",
            maintain_order="left",
        )
        .select([f"aux_{i}" for i in range(4)])
        .to_numpy()
        .astype(np.float32)
    )
    views = {
        "morgan": xf,
        "logfc": xn,
        "class": xc,
        "mixed": xm,
        "tree": xt,
        "class_logfc": np.column_stack([xc, xn]),
        "mixed_logfc": np.column_stack([xm, xn]),
        "tree_logfc": np.column_stack([xt, xn]),
    }
    individual = {}
    for candidate in selected:
        name = candidate["name"]
        booster = lgb.Booster(model_file=str(inp / record["models"][name]))
        individual[name] = booster.predict(views[candidate["view"]][ix])
        candidate_specs[name] = {"view": candidate["view"], "params": candidate["params"]}
    saved_predictions[fold] = individual
    reconstructed = sum(c["weight"] * individual[c["name"]] for c in selected)
    expected = fold_rows["prediction"].to_numpy()
    reconstruction_rows.append(
        {
            "fold": fold,
            "selected_models": len(selected),
            "best_step": selection["best_step"],
            "max_abs_error": float(np.max(np.abs(reconstructed - expected))),
            "mean_abs_error": float(np.mean(np.abs(reconstructed - expected))),
        }
    )
assert max(row["max_abs_error"] for row in reconstruction_rows) < 1e-10


def normalized_top(selected: list[dict], k: int, equal: bool) -> list[tuple[str, float]]:
    retained = selected[: min(k, len(selected))]
    if equal:
        return [(item["name"], 1 / len(retained)) for item in retained]
    total = sum(item["weight"] for item in retained)
    return [(item["name"], item["weight"] / total) for item in retained]


predictions = []
fold_score_rows = []
score_rows = []
max_k = cfg["max_models"]

# Direct truncation of each fold's saved inner-selected ensemble.
for weighting in ["weighted", "equal"]:
    for k in range(1, max_k + 1):
        prediction = np.empty(base.height)
        fit_count = 0
        for fold in range(5):
            mask = base["fold"].to_numpy() == fold
            chosen = normalized_top(selection_by_fold[fold], k, weighting == "equal")
            fit_count += len(chosen)
            prediction[mask] = sum(
                weight * saved_predictions[fold][name] for name, weight in chosen
            )
            fold_frame = base.filter(pl.col("fold") == fold)
            fold_ensemble = (prediction[mask] + other_sum[mask]) / len(ensemble_frames)
            fold_score_rows.append(
                {
                    "strategy": "nested_topk",
                    "weighting": weighting,
                    "k": k,
                    "fold": fold,
                    "retained_models": len(chosen),
                    "caruana_ST_RAE": score(fold_frame, prediction[mask]),
                    "ensemble_ST_RAE": score(fold_frame, fold_ensemble),
                }
            )
        predictions.append(("nested_topk", weighting, k, prediction, fit_count))

# Simulate the full-data recipe without using the held-out fold's selected weights.
for outer in range(5):
    fold_rows = base.filter(pl.col("fold") == outer)
    fold_target = target.filter((pl.col("fold") == outer) & ~pl.col("is_anchor"))
    outer_valid = (target["fold"].to_numpy() == outer) & ~target["is_anchor"].to_numpy()
    outer_train = ~outer_valid
    ix_train = target_indices[outer_train]
    ix_valid = np.array([index[c] for c in fold_rows["compound_id"]])
    y_train = target["value"].to_numpy()[outer_train]
    nested = pl.read_parquet(inp / fold_records[outer]["nested"])
    xn = (
        ids.join(
            nested.filter(pl.col("outer_fold") == outer),
            on="compound_id",
            validate="1:1",
            maintain_order="left",
        )
        .select([f"aux_{i}" for i in range(4)])
        .to_numpy()
        .astype(np.float32)
    )
    views = {
        "morgan": xf,
        "logfc": xn,
        "class": xc,
        "mixed": xm,
        "tree": xt,
        "class_logfc": np.column_stack([xc, xn]),
        "mixed_logfc": np.column_stack([xm, xn]),
        "tree_logfc": np.column_stack([xt, xn]),
    }
    aggregate = {}
    for source_fold, selected in selection_by_fold.items():
        if source_fold == outer:
            continue
        for item in selected:
            aggregate[item["name"]] = aggregate.get(item["name"], 0.0) + item["weight"] / 4
    ranked = sorted(aggregate.items(), key=lambda item: (-item[1], item[0]))
    model_predictions = {}
    for name, _weight in ranked[:max_k]:
        spec_item = candidate_specs[name]
        params = dict(spec_item["params"], random_state=cfg["seed"] + outer)
        model = lgb.LGBMRegressor(**params)
        model.fit(views[spec_item["view"]][ix_train], y_train)
        model_predictions[name] = model.predict(views[spec_item["view"]][ix_valid])
    for name in fixed_atomic_rows:
        assert name in model_predictions
        fixed_atomic_rows[name].append(
            fold_rows.select(keys + truth).with_columns(
                pl.Series("prediction", model_predictions[name]),
                pl.lit(name).alias("method"),
            )
        )
    for weighting in ["weighted", "equal"]:
        for k in range(1, max_k + 1):
            retained = ranked[: min(k, len(ranked))]
            if weighting == "equal":
                weights = [(name, 1 / len(retained)) for name, _ in retained]
            else:
                total = sum(weight for _, weight in retained)
                weights = [(name, weight / total) for name, weight in retained]
            p = sum(weight * model_predictions[name] for name, weight in weights)
            predictions.append(("crossfold_topk", weighting, k, outer, p, len(retained)))

# Join the fold-wise crossfold predictions into one OOF vector.
joined = []
for weighting in ["weighted", "equal"]:
    for k in range(1, max_k + 1):
        prediction = np.empty(base.height)
        fit_count = 0
        for item in predictions:
            if len(item) == 6 and item[:3] == ("crossfold_topk", weighting, k):
                _, _, _, fold, p, count = item
                mask = base["fold"].to_numpy() == fold
                prediction[mask] = p
                fit_count += count
                fold_frame = base.filter(pl.col("fold") == fold)
                fold_score_rows.append(
                    {
                        "strategy": "crossfold_topk",
                        "weighting": weighting,
                        "k": k,
                        "fold": fold,
                        "retained_models": count,
                        "caruana_ST_RAE": score(fold_frame, p),
                        "ensemble_ST_RAE": score(
                            fold_frame, (p + other_sum[mask]) / len(ensemble_frames)
                        ),
                    }
                )
        joined.append(("crossfold_topk", weighting, k, prediction, fit_count))
predictions = [item for item in predictions if len(item) == 5] + joined

rng = np.random.default_rng(cfg["bootstrap_seed"])
for strategy, weighting, k, prediction, fit_count in predictions:
    ensemble = (prediction + other_sum) / len(ensemble_frames)
    delta = score(base, ensemble) - score(base, original_ensemble)
    bootstrap = []
    for _ in range(cfg["bootstrap_replicates"]):
        sample = rng.integers(0, base.height, base.height)
        sampled = base[sample]
        bootstrap.append(
            score(sampled, ensemble[sample]) - score(sampled, original_ensemble[sample])
        )
    score_rows.append(
        {
            "strategy": strategy,
            "weighting": weighting,
            "k": k,
            "fold_fit_count": fit_count,
            "caruana_ST_RAE": score(base, prediction),
            "ensemble_ST_RAE": score(base, ensemble),
            "delta_vs_original_ensemble": delta,
            "delta_ci_low": float(np.quantile(bootstrap, 0.025)),
            "delta_ci_high": float(np.quantile(bootstrap, 0.975)),
        }
    )

pl.DataFrame(reconstruction_rows).write_csv(out / "reconstruction.csv")
pl.DataFrame(fold_score_rows).sort("strategy", "weighting", "k", "fold").write_csv(
    out / "fold-scores.csv"
)
pl.DataFrame(score_rows).sort("strategy", "weighting", "k").write_csv(out / "scores.csv")
prediction_frames = []
for strategy, weighting, k, prediction, _fit_count in predictions:
    prediction_frames.append(
        base.select(keys + truth).with_columns(
            pl.Series("caruana_prediction", prediction),
            pl.Series("ensemble_prediction", (prediction + other_sum) / len(ensemble_frames)),
            pl.lit(strategy).alias("strategy"),
            pl.lit(weighting).alias("weighting"),
            pl.lit(k).alias("k"),
        )
    )
pl.concat(prediction_frames).write_parquet(out / "predictions.parquet")
for name, frames in fixed_atomic_rows.items():
    pl.concat(frames).sort(keys).write_parquet(out / f"atomic-{name}.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "endpoint": endpoint,
            "original_caruana_ST_RAE": score(base, original_caruana),
            "ensemble_components": ensemble_names,
            "original_ensemble_ST_RAE": score(base, original_ensemble),
            "original_fold_selected_models": [len(selection_by_fold[fold]) for fold in range(5)],
            "original_fold_fit_count": sum(len(v) for v in selection_by_fold.values()),
            "deployment_distinct_models": cfg["deployment_distinct_models"],
            "reconstruction_max_abs_error": max(
                row["max_abs_error"] for row in reconstruction_rows
            ),
            "methods": {
                "nested_topk": "Keep each fold's K highest inner-selected weights.",
                "crossfold_topk": (
                    "For each held-out fold, average weights from the other four folds, "
                    "retain the top K candidates, and refit them on that fold's training labels."
                ),
                "weighted": "Renormalize retained Caruana weights to sum to one.",
                "equal": "Average retained candidate predictions equally.",
            },
        },
        indent=2,
    )
)
