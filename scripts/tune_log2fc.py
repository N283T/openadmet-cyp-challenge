"""Recorded Rustuna search; outer folds 0 and 1 never enter tuning."""

import gc
import json
import os
import traceback
from pathlib import Path

import polars as pl
import rustuna
import torch

from cyp_challenge.log2fc_pretrain import main as train
from cyp_challenge.storage import write_json

out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
base = cfg["base"]
study = rustuna.create_study(
    direction="minimize", sampler=rustuna.samplers.TPESampler(seed=42, n_startup_trials=6)
)
records = []


def objective(trial):
    trial_out = out / f"trial-{trial.number:03d}"
    trial_out.mkdir()
    params = (
        {}
        if trial.number == 0
        else {
            "learning_rate": trial.suggest_float("head_lr", 2e-5, 1e-3, log=True),
            "lr_ratio": 1.0,
            "encoder_lr_scale": trial.suggest_float("encoder_lr_scale", 0.03, 1.0, log=True),
            "weight_decay": trial.suggest_float("weight_decay", 1e-6, 0.05, log=True),
            "ffn_dropout": trial.suggest_float("ffn_dropout", 0.0, 0.5),
            "batch_norm": trial.suggest_categorical("batch_norm", [True, False]),
            "scheduler": trial.suggest_categorical("scheduler", ["plateau", "noam"]),
            "optimizer": "adamw",
            "max_epochs": 80,
            "patience": 15,
            "lr_patience": 3,
        }
    )
    config = {**base, **params, "tuning_only": True}
    # For noam, sampled rate denotes peak; retain a 10x warmup ramp.
    if trial.number and config["scheduler"] == "noam":
        config["learning_rate"] /= 10
        config["lr_ratio"] = 10.0
    config_path = trial_out / "config.json"
    write_json(config_path, config)
    os.environ["CYP_OUTPUT_DIR"] = str(trial_out)
    os.environ["CYP_CONFIG_PATH"] = str(config_path)
    try:
        train()
        result = json.loads((trial_out / "summary.json").read_text())
        value = result["best_validation_mse_z"]
        records.append(
            {
                "trial": trial.number,
                "value": value,
                "best_epoch": result["best_epoch"],
                "epochs": result["epochs"],
                "params": params,
                "config": config,
            }
        )
        with (out / "trials.jsonl").open("a") as stream:
            stream.write(json.dumps(records[-1]) + "\n")
        print(f"RUSTUNA trial={trial.number} value={value}", flush=True)
        return value
    except Exception:
        (trial_out / "failure.txt").write_text(traceback.format_exc())
        raise
    finally:
        gc.collect()
        torch.cuda.empty_cache()


study.optimize(objective, n_trials=cfg["n_trials"])
write_json(out / "trials.json", records)
winner = min(records, key=lambda r: r["value"])
write_json(out / "selected.json", winner)
pl.DataFrame(
    [{k: r[k] for k in ["trial", "value", "best_epoch", "epochs"]} for r in records]
).write_csv(out / "trials.csv")
print("SELECTED", json.dumps(winner), flush=True)
