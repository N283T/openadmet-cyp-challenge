"""Fold-excluded recipe selection from the fixed robust-bound blend screen."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def st_rae(y: np.ndarray, lo: np.ndarray, hi: np.ndarray, pred: np.ndarray) -> float:
    center = y.mean()
    return float(
        (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
        / (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    )


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    all_oof = pl.read_parquet(inp / "blend-oof.parquet").filter(
        (pl.col("placement") == cfg["placement"]) & (pl.col("blend_weight") > 0)
    )
    keys = ["robust_delta", "constraint_weight", "blend_weight"]
    rows = []
    selected = []
    for split in ["parent", "scaffold"]:
        split_oof = all_oof.filter(pl.col("split_id") == split)
        for held_fold in range(5):
            train = split_oof.filter(pl.col("fold") != held_fold)
            scored = []
            for recipe, part in train.group_by(keys):
                if part.height < 1000:
                    raise ValueError("Insufficient recipe coverage")
                scored.append(
                    (
                        st_rae(
                            *[
                                part[col].to_numpy()
                                for col in ["y_true", "conf_low", "conf_high", "prediction"]
                            ]
                        ),
                        recipe,
                    )
                )
            train_score, recipe = min(scored, key=lambda item: (item[0], item[1]))
            held = split_oof.filter(
                (pl.col("fold") == held_fold)
                & (pl.col("robust_delta") == recipe[0])
                & (pl.col("constraint_weight") == recipe[1])
                & (pl.col("blend_weight") == recipe[2])
            )
            if held.height < 200:
                raise ValueError("Incomplete held-out fold")
            selected.append(held)
            rows.append(
                {
                    "split_id": split,
                    "held_fold": held_fold,
                    "robust_delta": recipe[0],
                    "constraint_weight": recipe[1],
                    "blend_weight": recipe[2],
                    "train_st_rae": train_score,
                    "held_st_rae": st_rae(
                        *[
                            held[col].to_numpy()
                            for col in ["y_true", "conf_low", "conf_high", "prediction"]
                        ]
                    ),
                }
            )
    oof = pl.concat(selected).sort("split_id", "compound_id")
    if oof.select("split_id", "compound_id").is_duplicated().any():
        raise ValueError("Duplicate selected OOF key")
    summary = []
    for (split,), part in oof.group_by("split_id"):
        summary.append(
            {
                "split_id": split,
                "n": part.height,
                "ST_RAE": st_rae(
                    *[
                        part[col].to_numpy()
                        for col in ["y_true", "conf_low", "conf_high", "prediction"]
                    ]
                ),
            }
        )
    pl.DataFrame(rows).sort("split_id", "held_fold").write_csv(out / "selection.csv")
    pl.DataFrame(summary).sort("split_id").write_csv(out / "metrics.csv")
    oof.write_parquet(out / "oof.parquet")
    (out / "summary.json").write_text(json.dumps({"config": cfg}, indent=2) + "\n")


if __name__ == "__main__":
    main()
