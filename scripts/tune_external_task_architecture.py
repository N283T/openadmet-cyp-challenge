"""Architecture screen for the fixed intrinsic CYP continuous task."""

import gc
import json
import os
import shutil
from pathlib import Path

import torch
import train_external_medium as medium

out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg_path = Path(os.environ["CYP_CONFIG_PATH"])
cfg = json.loads(cfg_path.read_text())
medium.trainer.AssayModel = medium.MediumAssayModel

results = []
for index, variant in enumerate(cfg["variants"]):
    trial_dir = out / f"trial-{index:02d}-{variant['name']}"
    trial_dir.mkdir()
    trial_cfg = {**cfg["base"], **{k: v for k, v in variant.items() if k != "name"}}
    trial_cfg["variant"] = variant["name"]
    trial_config_path = trial_dir / "config.json"
    trial_config_path.write_text(json.dumps(trial_cfg, indent=2) + "\n")
    os.environ["CYP_OUTPUT_DIR"] = str(trial_dir)
    os.environ["CYP_CONFIG_PATH"] = str(trial_config_path)
    medium.trainer.main()
    summary = json.loads((trial_dir / "summary.json").read_text())
    result = {
        "variant": variant["name"],
        "validation_loss": summary["best_validation_loss"],
        "best_epoch": summary["best_epoch"],
        "epochs": summary["epochs"],
        "assays": summary["assays"],
        "output_heads": summary["output_heads"],
        "tasks": summary["tasks"],
        "parameters": variant,
    }
    results.append(result)
    (out / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(result), flush=True)
    gc.collect()
    torch.cuda.empty_cache()

os.environ["CYP_OUTPUT_DIR"] = str(out)
os.environ["CYP_CONFIG_PATH"] = str(cfg_path)
winner = min(results, key=lambda row: row["validation_loss"])
winner_dir = next(path for path in out.glob("trial-*") if path.name.endswith(winner["variant"]))
shutil.copy2(winner_dir / "best.pt", out / "selected-best.pt")
(out / "selection.json").write_text(
    json.dumps({"protocol": cfg["protocol"], "winner": winner}, indent=2) + "\n"
)
print(json.dumps({"winner": winner}, indent=2), flush=True)
