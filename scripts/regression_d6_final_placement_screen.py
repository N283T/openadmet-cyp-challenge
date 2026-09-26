"""Test label-free final placement of the selected D6 blend on both OOF splits."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata, spearmanr


def st_rae(frame: pl.DataFrame) -> float:
    y = frame["y_true"].to_numpy()
    p = frame["prediction"].to_numpy()
    lo = frame["conf_low"].to_numpy()
    hi = frame["conf_high"].to_numpy()
    center = y.mean()
    return float(
        (np.maximum(lo - p, 0) + np.maximum(p - hi, 0)).sum()
        / (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    )


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    selected = pl.read_parquet(inp / "blend-oof.parquet").filter(
        (pl.col("endpoint") == cfg["endpoint"])
        & (pl.col("robust_delta") == cfg["robust_delta"])
        & (pl.col("constraint_weight") == cfg["constraint_weight"])
        & (pl.col("placement") == cfg["placement"])
        & pl.col("blend_weight").is_in([0.0, cfg["blend_weight"]])
    )
    outputs = []
    moments = []
    for split in ["parent", "scaffold"]:
        for fold in range(5):
            part = selected.filter((pl.col("split_id") == split) & (pl.col("fold") == fold))
            base = part.filter(pl.col("blend_weight") == 0.0).sort("compound_id")
            candidate = part.filter(pl.col("blend_weight") == cfg["blend_weight"]).sort(
                "compound_id"
            )
            if base.height < 200 or candidate.height != base.height:
                raise ValueError("Incomplete fold pair")
            for key in ["compound_id", "endpoint", "y_true", "conf_low", "conf_high"]:
                if not base[key].equals(candidate[key]):
                    raise ValueError(f"Mismatched {split}/{fold}/{key}")
            x = candidate["prediction"].to_numpy()
            reference = base["prediction"].to_numpy()
            ranks = (rankdata(x, method="average") - 0.5) / len(x)
            methods = {
                "identity": x,
                "restore_control_sd": (x - x.mean()) / x.std() * reference.std() + reference.mean(),
                "match_control_quantiles": np.quantile(reference, ranks),
            }
            for method in cfg["methods"]:
                pred = methods[method]
                outputs.append(
                    candidate.select(
                        "compound_id",
                        "endpoint",
                        "split_id",
                        "fold",
                        "y_true",
                        "conf_low",
                        "conf_high",
                    ).with_columns(pl.lit(method).alias("method"), pl.Series("prediction", pred))
                )
                moments.append(
                    {
                        "split_id": split,
                        "fold": fold,
                        "method": method,
                        "mean": float(pred.mean()),
                        "sd": float(pred.std()),
                        "control_mean": float(reference.mean()),
                        "control_sd": float(reference.std()),
                    }
                )
    oof = pl.concat(outputs).sort("split_id", "method", "compound_id")
    if oof.select("split_id", "method", "compound_id").is_duplicated().any():
        raise ValueError("Duplicate OOF key")
    metrics = []
    for (split, method), part in oof.group_by("split_id", "method"):
        metrics.append(
            {
                "split_id": split,
                "method": method,
                "n": part.height,
                "ST_RAE": st_rae(part),
                "Spearman": float(spearmanr(part["y_true"], part["prediction"]).statistic),
            }
        )
    oof.write_parquet(out / "oof.parquet")
    pl.DataFrame(metrics).sort("split_id", "ST_RAE").write_csv(out / "metrics.csv")
    pl.DataFrame(moments).sort("split_id", "fold", "method").write_csv(out / "moments.csv")
    (out / "summary.json").write_text(json.dumps({"config": cfg}, indent=2) + "\n")


if __name__ == "__main__":
    main()
