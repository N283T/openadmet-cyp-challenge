"""Rustuna search for one project-owned encoder checkpoint and raw embedding."""

import gc
import importlib.util
import json
import os
import shutil
import traceback
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
import rustuna
import torch
import train_external_medium as medium

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
medium.trainer.AssayModel = medium.MediumAssayModel


def input_path(name: str) -> Path:
    exact = inp / name
    paths = [exact] if exact.exists() else [
        path for path in inp.iterdir() if path.name.partition("-")[2] == name
    ]
    if len(paths) != 1:
        raise ValueError((name, paths))
    return paths[0]


spec = importlib.util.spec_from_file_location(
    "official", input_path("custom_scoring_functions.py")
)
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
assignments = pl.read_parquet(input_path("parent-assignments.parquet"))
labels = pl.read_parquet(input_path("measurements.parquet")).filter(
    pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
)
jeremy = pl.read_parquet(input_path("jeremy-oof.parquet"))
endpoints = sorted(labels["endpoint"].unique())
assert len(endpoints) == 4


def score_embedding(path: Path, validation_fold: int) -> tuple[float, list[dict], pl.DataFrame]:
    emb = pl.read_parquet(path)
    columns = sorted(column for column in emb.columns if column.startswith("emb_"))
    assert columns and emb["compound_id"].n_unique() == emb.height
    metrics, rows = [], []
    for endpoint in endpoints:
        table = (
            labels.filter(pl.col("endpoint") == endpoint)
            .join(assignments, on="compound_id", validate="1:1")
            .join(emb, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        x = table.select(columns).to_numpy()
        y = table["value"].to_numpy()
        valid = (table["fold"].to_numpy() == validation_fold) & ~table[
            "is_anchor"
        ].to_numpy()
        model = lgb.LGBMRegressor(
            objective="regression_l1",
            n_estimators=300,
            num_leaves=15,
            min_child_samples=20,
            learning_rate=0.05,
            reg_lambda=1.0,
            n_jobs=4,
            random_state=42 + validation_fold,
            verbosity=-1,
            deterministic=True,
            force_col_wise=True,
        )
        model.fit(x[~valid], y[~valid])
        prediction = model.predict(x[valid])
        selected = table.filter(pl.Series(valid))
        score = float(
            scorer.rae_soft_threshold_absolute_error(
                selected["value"].to_numpy(),
                prediction,
                y_true_lower=selected["conf_low"].to_numpy(),
                y_true_upper=selected["conf_high"].to_numpy(),
            )
        )
        metrics.append({"endpoint": endpoint, "fold": validation_fold, "ST_RAE": score})
        rows.append(
            selected.select(
                "compound_id",
                "endpoint",
                "fold",
                "is_neighbor",
                "conf_low",
                "conf_high",
                pl.col("value").alias("y_true"),
            ).with_columns(pl.Series("prediction", prediction))
        )
    macro = float(np.mean([row["ST_RAE"] for row in metrics]))
    return macro, metrics, pl.concat(rows)


def reference_scores(fold: int) -> tuple[float, list[dict]]:
    metrics = []
    for endpoint in endpoints:
        table = jeremy.filter(
            (pl.col("endpoint") == endpoint) & (pl.col("fold") == fold)
        )
        value = float(
            scorer.rae_soft_threshold_absolute_error(
                table["y_true"].to_numpy(),
                table["prediction"].to_numpy(),
                y_true_lower=table["conf_low"].to_numpy(),
                y_true_upper=table["conf_high"].to_numpy(),
            )
        )
        metrics.append({"endpoint": endpoint, "fold": fold, "ST_RAE": value})
    return float(np.mean([row["ST_RAE"] for row in metrics])), metrics


control = {
    "width": 600,
    "depth": 3,
    "aggregation": "norm",
    "aggregation_norm": 100.0,
    "ffn_hidden_layers": 3,
    "dropout": 0.1,
    "learning_rate": 1e-3,
    "min_learning_rate": 1e-4,
    "weight_decay": 1e-4,
    "batch_size": 128,
    "rank_min_gap": 0.25,
    "continuous_weight": 1.0,
    "ranking_weight": 1.0,
    "classification_weight": 1.0,
    "scaling_pseudocount": 20.0,
    "scaling_minimum_std": 0.3,
}
records = []
study = rustuna.create_study(
    direction="minimize",
    sampler=rustuna.samplers.TPESampler(
        seed=cfg["sampler_seed"], n_startup_trials=6
    ),
)


def objective(trial):
    trial_dir = out / f"trial-{trial.number:03d}"
    trial_dir.mkdir()
    if trial.number == 0:
        params = control
    else:
        learning_rate = trial.suggest_float("learning_rate", 2e-4, 2e-3, log=True)
        params = {
            "width": trial.suggest_categorical(
                "width", [448, 512, 600, 768, 896, 1024]
            ),
            "depth": trial.suggest_int("depth", 2, 5),
            "aggregation": "norm",
            "aggregation_norm": trial.suggest_float(
                "aggregation_norm", 30.0, 300.0, log=True
            ),
            "ffn_hidden_layers": trial.suggest_int("ffn_hidden_layers", 1, 4),
            "dropout": trial.suggest_float("dropout", 0.0, 0.3),
            "learning_rate": learning_rate,
            "min_learning_rate": learning_rate
            * trial.suggest_float("min_lr_ratio", 0.03, 0.3, log=True),
            "weight_decay": trial.suggest_float("weight_decay", 1e-6, 1e-2, log=True),
            "batch_size": trial.suggest_categorical("batch_size", [64, 128, 256]),
            "rank_min_gap": trial.suggest_float("rank_min_gap", 0.1, 0.6),
            "continuous_weight": trial.suggest_float(
                "continuous_weight", 0.4, 2.5, log=True
            ),
            "ranking_weight": trial.suggest_float("ranking_weight", 0.2, 2.5, log=True),
            "classification_weight": trial.suggest_float(
                "classification_weight", 0.2, 2.5, log=True
            ),
            "scaling_pseudocount": trial.suggest_float(
                "scaling_pseudocount", 5.0, 60.0, log=True
            ),
            "scaling_minimum_std": trial.suggest_float(
                "scaling_minimum_std", 0.15, 0.7
            ),
        }
    trial_cfg = {**cfg["base"], **params, "trial": trial.number}
    trial_config_path = trial_dir / "config.json"
    trial_config_path.write_text(json.dumps(trial_cfg, indent=2) + "\n")
    os.environ["CYP_OUTPUT_DIR"] = str(trial_dir)
    os.environ["CYP_CONFIG_PATH"] = str(trial_config_path)
    try:
        medium.trainer.main()
        macro, metrics, predictions = score_embedding(
            trial_dir / "external-embeddings.parquet", cfg["tuning_fold"]
        )
        predictions.write_parquet(trial_dir / "tuning-fold-oof.parquet")
        training = json.loads((trial_dir / "summary.json").read_text())
        record = {
            "trial": trial.number,
            "score": macro,
            "metrics": metrics,
            "params": params,
            "training": {
                key: training[key]
                for key in [
                    "assays",
                    "tasks",
                    "best_epoch",
                    "epochs",
                    "best_validation_loss",
                    "embedding_dimensions",
                ]
            },
            "status": "complete",
        }
    except Exception:
        (trial_dir / "failure.txt").write_text(traceback.format_exc())
        record = {
            "trial": trial.number,
            "score": 10.0,
            "metrics": [],
            "params": params,
            "status": "failed",
        }
    finally:
        os.environ["CYP_OUTPUT_DIR"] = str(out)
        os.environ["CYP_CONFIG_PATH"] = str(Path(cfg["config_path"]))
        gc.collect()
        torch.cuda.empty_cache()
    records.append(record)
    with (out / "trials.jsonl").open("a") as stream:
        stream.write(json.dumps(record) + "\n")
    print(
        f"RUSTUNA trial={trial.number} score={record['score']} status={record['status']}",
        flush=True,
    )
    return record["score"]


cfg["config_path"] = os.environ["CYP_CONFIG_PATH"]
study.optimize(objective, n_trials=cfg["n_trials"])
completed = [record for record in records if record["status"] == "complete"]
if not completed:
    raise RuntimeError("Every Rustuna trial failed")
winner = min(completed, key=lambda record: record["score"])
winner_dir = out / f"trial-{winner['trial']:03d}"
for source, destination in [
    ("best.pt", "selected-best.pt"),
    ("external-embeddings.parquet", "selected-embeddings.parquet"),
    ("architecture.json", "selected-architecture.json"),
    ("summary.json", "selected-training-summary.json"),
    ("config.json", "selected-training-config.json"),
]:
    shutil.copy2(winner_dir / source, out / destination)
holdout_macro, holdout_metrics, holdout_predictions = score_embedding(
    out / "selected-embeddings.parquet", cfg["untouched_holdout_fold"]
)
holdout_predictions.write_parquet(out / "holdout-oof.parquet")
jeremy_holdout_macro, jeremy_holdout_metrics = reference_scores(
    cfg["untouched_holdout_fold"]
)
selection = {
    "winner": winner,
    "holdout_macro_ST_RAE": holdout_macro,
    "holdout_metrics": holdout_metrics,
    "jeremy_holdout_macro_ST_RAE": jeremy_holdout_macro,
    "jeremy_holdout_metrics": jeremy_holdout_metrics,
    "holdout_delta": holdout_macro - jeremy_holdout_macro,
    "holdout_was_unread_during_selection": True,
}
(out / "selected.json").write_text(json.dumps(selection, indent=2) + "\n")
(out / "trials.json").write_text(json.dumps(records, indent=2) + "\n")
pl.DataFrame(
    [
        {
            "trial": record["trial"],
            "score": record["score"],
            "status": record["status"],
        }
        for record in records
    ]
).sort("score").write_csv(out / "trials.csv")
print("SELECTED", json.dumps(selection), flush=True)
