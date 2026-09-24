"""Freeze exact accepted classification vectors and endpoint LB transitions."""

import json
import os
from datetime import datetime
from pathlib import Path

import polars as pl

inputs = Path(os.environ["CYP_INPUT_DIR"])
outputs = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoints = cfg["endpoint_columns"]
states = cfg["states"]
if endpoints != ["CYP2D6_is_TDI", "CYP3A4_is_TDI"] or len(states) != 7:
    raise ValueError("Unexpected frozen historical state definition")

base_ids = None
base_smiles = None
state_rows = []
transition_rows = []
swap_rows = []
previous = None
for state in states:
    path = inputs / state["submission"]
    package = pl.read_parquet(path).sort("Molecule_Name")
    if package.height != 750 or package["Molecule_Name"].n_unique() != 750:
        raise ValueError(f"Invalid submission vector: {state['date']}")
    ids = package["Molecule_Name"].to_list()
    smiles = package["SMILES"].to_list()
    if base_ids is None:
        base_ids, base_smiles = ids, smiles
    elif ids != base_ids or smiles != base_smiles:
        raise ValueError(f"Submission identity mismatch: {state['date']}")
    calls = {}
    scores = {}
    for endpoint, key in zip(endpoints, ["d6_lb", "a4_lb"], strict=True):
        calls[endpoint] = package[endpoint].cast(pl.Boolean).to_list()
        leaderboard = pl.read_csv(inputs / state[key]).filter(
            pl.col("Username") == cfg["username"]
        )
        if leaderboard.height != 1:
            raise ValueError(f"Missing or duplicate own LB row: {state['date']} {endpoint}")
        row = leaderboard.row(0, named=True)
        submitted = datetime.strptime(row["Submitted"], "%Y-%m-%d %H:%M UTC")
        if submitted.date().isoformat() != state["date"]:
            raise ValueError(f"LB timestamp mismatch: {state['date']} {endpoint}")
        scores[endpoint] = float(row["MCC"])
        state_rows.append(
            {
                "date": state["date"],
                "endpoint": endpoint,
                "submission_alias": state["submission"],
                "leaderboard_alias": state[key],
                "submitted_utc": row["Submitted"],
                "mcc": scores[endpoint],
                "mcc_std": float(row["MCC std"]),
                "positive_calls": sum(calls[endpoint]),
            }
        )
    if previous is not None:
        for endpoint in endpoints:
            old, new = previous["calls"][endpoint], calls[endpoint]
            removed = [id_ for id_, a, b in zip(ids, old, new, strict=True) if a and not b]
            added = [id_ for id_, a, b in zip(ids, old, new, strict=True) if b and not a]
            swap_rows.extend(
                {"from_date": previous["date"], "to_date": state["date"],
                 "endpoint": endpoint, "molecule_name": id_, "change": change}
                for change, group in (("removed", removed), ("added", added))
                for id_ in group
            )
            delta = scores[endpoint] - previous["scores"][endpoint]
            transition_rows.append(
                {
                    "from_date": previous["date"],
                    "to_date": state["date"],
                    "endpoint": endpoint,
                    "lb_mcc_delta": delta,
                    "direction": "gain" if delta > 1e-9 else "loss" if delta < -1e-9 else "tie",
                    "old_calls": sum(old),
                    "new_calls": sum(new),
                    "removed_count": len(removed),
                    "added_count": len(added),
                }
            )
    previous = {"date": state["date"], "calls": calls, "scores": scores}

pl.DataFrame(state_rows).write_csv(outputs / "states.csv")
pl.DataFrame(transition_rows).write_csv(outputs / "transitions.csv")
pl.DataFrame(swap_rows).write_parquet(outputs / "changed-calls.parquet")
(outputs / "summary.json").write_text(
    json.dumps(
        {
            "distinct_states": len(states),
            "endpoint_transitions": len(transition_rows),
            "molecule_count": len(base_ids),
            "d6_directions": [
                row["direction"]
                for row in transition_rows
                if row["endpoint"] == "CYP2D6_is_TDI"
            ],
            "a4_directions": [
                row["direction"]
                for row in transition_rows
                if row["endpoint"] == "CYP3A4_is_TDI"
            ],
        },
        indent=2,
    )
    + "\n"
)
