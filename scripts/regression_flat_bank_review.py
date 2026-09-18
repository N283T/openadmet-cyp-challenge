"""Run auditable target-specific variants of the flat regression banks."""

import json
import os
import subprocess
from pathlib import Path

import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
scratch = Path(os.environ["TMPDIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
bank_script = Path(__file__).with_name("regression_flat_ensemble_bank.py")

summary_rows = []
for target, pool in cfg["pools"].items():
    available = {name: file for name, file in pool["components"]}
    for variant, names in pool["variants"].items():
        if len(names) != len(set(names)) or not set(names).issubset(available):
            raise ValueError(f"Invalid members for {target}/{variant}")
        variant_cfg = {
            "endpoints": [pool["endpoint"]],
            "components": [
                {"name": name, "file": available[name]} for name in names
            ],
            "current_control": cfg["current_control"],
            "scorer": cfg["scorer"],
            "ridge_alpha": cfg["ridge_alpha"],
            "simplex_lambda": cfg["simplex_lambda"],
            "caruana_steps": cfg["caruana_steps"],
            "protocol": f"Flat-bank review {target}/{variant}",
        }
        config_path = scratch / f"{target}-{variant}.json"
        config_path.write_text(json.dumps(variant_cfg, indent=2) + "\n")
        variant_out = out / target / variant
        variant_out.mkdir(parents=True)
        env = dict(os.environ)
        env["CYP_CONFIG_PATH"] = str(config_path)
        env["CYP_OUTPUT_DIR"] = str(variant_out)
        subprocess.run(["python", str(bank_script)], check=True, env=env)
        metrics = pl.read_csv(variant_out / "metrics.csv")
        for row in metrics.iter_rows(named=True):
            summary_rows.append(
                {
                    "target": target,
                    "variant": variant,
                    "members": len(names),
                    **row,
                }
            )

summary = pl.DataFrame(summary_rows).sort("target", "ST_RAE")
summary.write_csv(out / "variant-metrics.csv")
best = summary.group_by("target", maintain_order=True).head(1)
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "best_by_target": best.to_dicts(),
            "selection_note": (
                "Exploratory variants chosen after inspecting the same OOF predictions; "
                "use fold stability and simplicity, not the minimum alone."
            ),
        },
        indent=2,
    )
    + "\n"
)
print(best, flush=True)
