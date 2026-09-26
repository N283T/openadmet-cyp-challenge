"""Paired group-bootstrap audit of one fixed D6 primary-bound blend."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr


def score(y: np.ndarray, lo: np.ndarray, hi: np.ndarray, pred: np.ndarray) -> float:
    center = y.mean()
    numerator = (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
    denominator = (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    return float(numerator / denominator)


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    blend = pl.read_parquet(inp / "blend-oof.parquet").filter(
        (pl.col("endpoint") == cfg["endpoint"])
        & (pl.col("constraint_weight") == cfg["constraint_weight"])
        & (pl.col("placement") == cfg["placement"])
        & pl.col("blend_weight").is_in([0.0, cfg["candidate_blend_weight"]])
    )
    rows = []
    draws = []
    influence = []
    for split in ["parent", "scaffold"]:
        assignment = pl.read_parquet(inp / f"{split}-assignments.parquet").select(
            "compound_id", "fold", "group", "is_anchor"
        )
        frame = blend.filter(pl.col("split_id") == split).join(
            assignment, on=["compound_id", "fold"], validate="m:1"
        )
        if frame.height != 2 * 1457 or frame["is_anchor"].any():
            raise ValueError(f"Invalid {split} assignment coverage")
        base = frame.filter(pl.col("blend_weight") == 0).sort("compound_id")
        candidate = frame.filter(pl.col("blend_weight") == cfg["candidate_blend_weight"]).sort(
            "compound_id"
        )
        if base.height != 1457 or candidate.height != 1457:
            raise ValueError("Missing OOF arm")
        for key in ["compound_id", "endpoint", "fold", "group", "y_true", "conf_low", "conf_high"]:
            if not base[key].equals(candidate[key]):
                raise ValueError(f"Mismatched {split}/{key}")
        y = base["y_true"].to_numpy()
        lo = base["conf_low"].to_numpy()
        hi = base["conf_high"].to_numpy()
        p0 = base["prediction"].to_numpy()
        p1 = candidate["prediction"].to_numpy()
        groups = base["group"].to_numpy()
        folds = base["fold"].to_numpy()
        for fold in [-1, 0, 1, 2, 3, 4]:
            mask = np.ones(len(y), dtype=bool) if fold == -1 else folds == fold
            s0 = score(y[mask], lo[mask], hi[mask], p0[mask])
            s1 = score(y[mask], lo[mask], hi[mask], p1[mask])
            rows.append(
                {
                    "split_id": split,
                    "fold": fold,
                    "n": int(mask.sum()),
                    "base_st_rae": s0,
                    "candidate_st_rae": s1,
                    "delta_st_rae": s1 - s0,
                    "base_spearman": float(spearmanr(y[mask], p0[mask]).statistic),
                    "candidate_spearman": float(spearmanr(y[mask], p1[mask]).statistic),
                }
            )
        unique, inverse = np.unique(groups, return_inverse=True)
        members = [np.flatnonzero(inverse == index) for index in range(len(unique))]
        rng = np.random.default_rng(cfg["seed"] + (0 if split == "parent" else 1))
        for draw in range(cfg["bootstrap_draws"]):
            sample_groups = rng.integers(0, len(unique), size=len(unique))
            sample = np.concatenate([members[index] for index in sample_groups])
            draws.append(
                {
                    "split_id": split,
                    "draw": draw,
                    "delta_st_rae": score(y[sample], lo[sample], hi[sample], p1[sample])
                    - score(y[sample], lo[sample], hi[sample], p0[sample]),
                }
            )
        soft0 = np.maximum(lo - p0, 0) + np.maximum(p0 - hi, 0)
        soft1 = np.maximum(lo - p1, 0) + np.maximum(p1 - hi, 0)
        for index, group in enumerate(unique):
            selected = members[index]
            influence.append(
                {
                    "split_id": split,
                    "group": str(group),
                    "n": len(selected),
                    "delta_soft_error_sum": float((soft1[selected] - soft0[selected]).sum()),
                }
            )
    boot = pl.DataFrame(draws)
    summary = []
    for split, part in boot.group_by("split_id"):
        vals = part["delta_st_rae"].to_numpy()
        summary.append(
            {
                "split_id": split[0],
                "draws": len(vals),
                "q025": float(np.quantile(vals, 0.025)),
                "median": float(np.median(vals)),
                "q975": float(np.quantile(vals, 0.975)),
                "p_improve": float(np.mean(vals < 0)),
            }
        )
    pl.DataFrame(rows).sort("split_id", "fold").write_csv(out / "fold-metrics.csv")
    boot.write_parquet(out / "bootstrap.parquet")
    pl.DataFrame(summary).sort("split_id").write_csv(out / "bootstrap-summary.csv")
    pl.DataFrame(influence).sort("split_id", "delta_soft_error_sum").write_csv(
        out / "group-influence.csv"
    )
    (out / "summary.json").write_text(json.dumps({"config": cfg}, indent=2) + "\n")


if __name__ == "__main__":
    main()
