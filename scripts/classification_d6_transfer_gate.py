"""Apply a conservative complete-pool, dual-split D6 replacement evidence gate."""

import json
import os
from pathlib import Path

import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
checks = []
for strategy in cfg["strategies"]:
    for mode in ("fixed", "nested"):
        table = pl.read_csv(inp / cfg[f"{strategy}_{mode}"])
        selected = table.filter(
            (pl.col("strategy") == strategy)
            & (pl.col("scope") == cfg["scope"])
            & ((pl.col("prior") - cfg["target_prior"]).abs() < 1e-10)
        )
        if mode == "fixed":
            selected = selected.filter((pl.col("rate") - cfg["historical_rate"]).abs() < 1e-10)
        if selected.height != 1:
            raise ValueError(f"Expected one {strategy} {mode} full-pool result")
        row = selected.row(0, named=True)
        if row["pool_rows"] != cfg["expected_pool_rows"]:
            raise ValueError(f"Incomplete molecular pool for {strategy} {mode}")
        lower = float(row["delta_ci025"])
        checks.append(
            {
                "strategy": strategy,
                "decision": mode,
                "candidate_minus_control": float(row["candidate_minus_control"]),
                "lower_ci": lower,
                "upper_ci": float(row["delta_ci975"]),
                "positive_lower_interval": lower > cfg["minimum_lower_ci"],
            }
        )
status = (
    "ELIGIBLE_FOR_PROBE"
    if all(check["positive_lower_interval"] for check in checks)
    else "INSUFFICIENT_EVIDENCE"
)
pl.DataFrame(checks).write_csv(out / "checks.csv")
(out / "decision.json").write_text(
    json.dumps(
        {
            "status": status,
            "checks": checks,
            "minimum_lower_ci": cfg["minimum_lower_ci"],
            "meaning": (
                "conservative evidence requirement for considering a new blind probe; "
                "not a guarantee of leaderboard gain"
            ),
        },
        indent=2,
    )
    + "\n"
)
print(status, flush=True)
