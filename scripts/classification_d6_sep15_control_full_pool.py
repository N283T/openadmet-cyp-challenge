"""Recompose the September 15 D6 control from complete-pool OOF atoms."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
sources = {
    key: pl.read_parquet(inp / cfg[key]) for key in ["lgbm", "pfn_pca", "pfn_raw", "graph", "fixed"]
}
members = [
    ("lgbm", "morgan_plain"),
    ("lgbm", "mordred_plain"),
    ("lgbm", "fp_mordred_plain"),
    ("lgbm", "external_class_plain"),
    ("lgbm", "external_mixed_plain"),
    ("pfn_pca", "external_mixed_pca256"),
    ("pfn_raw", "external_mixed_raw600"),
    ("graph", "tdi_only"),
    ("fixed", "logistic_c0.05"),
    ("fixed", "rbf_svm_c1"),
    ("fixed", "complement_nb_a1"),
]
parts = []
for source, model in members:
    part = (
        sources[source]
        .filter(
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("split_id") == cfg["split_id"])
            & (pl.col("model") == model)
        )
        .sort("compound_id")
    )
    if part.height != 6146 or part["compound_id"].n_unique() != 6146:
        raise ValueError(f"Incomplete control atom: {source} {model}")
    parts.append(part)
ref = parts[0]
for part in parts[1:]:
    if not (
        part["compound_id"].equals(ref["compound_id"])
        and part["fold"].equals(ref["fold"])
        and part["y_true"].equals(ref["y_true"])
    ):
        raise ValueError("Control atom identity, fold or label disagreement")
folds = ref["fold"].to_numpy()
raw = np.column_stack([part["score"].to_numpy() for part in parts])
if not np.isfinite(raw).all():
    raise ValueError("Nonfinite control atom score")


def fold_rank(values: np.ndarray) -> np.ndarray:
    ranked = np.empty(len(values))
    for fold in np.unique(folds):
        mask = folds == fold
        ranked[mask] = rankdata(values[mask], method="average") / (mask.sum() + 1)
    return ranked


ranked = np.column_stack([fold_rank(raw[:, column]) for column in range(raw.shape[1])])
scores = {}
for mode in ("raw", "rank"):
    atoms = raw if mode == "raw" else ranked
    mixed = (atoms[:, 5] + atoms[:, 6]) / 2
    core = fold_rank(np.mean(np.column_stack([atoms[:, :5], mixed, atoms[:, 7]]), axis=1))
    logistic = ranked[:, 8]
    rbf = ranked[:, 9]
    complement = ranked[:, 10]
    live = (2 * core + logistic) / 3
    score = ((2 * live + rbf) / 3 + (2 * live + complement) / 3) / 2
    scores[f"sep15_control_{mode}_oof_analog"] = score
result = pl.concat(
    [
        ref.select("compound_id", "fold", "group", "y_true", "endpoint", "split_id").with_columns(
            pl.lit(name).alias("model"), pl.Series("score", score)
        )
        for name, score in scores.items()
    ]
)
result.write_parquet(out / "all-pool-oof.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "models": cfg["models"],
            "raw_analog": (
                "seven core raw scores, matching the blind package composition "
                "before rank; tails fold-ranked"
            ),
            "rank_analog": (
                "rank each core atom within the held fold before averaging, "
                "matching historical labelled OOF development"
            ),
            "scope": (
                "OOF analogs of September 15 fixed full-data recipe; "
                "the September 15 blind vector itself is unchanged"
            ),
        },
        indent=2,
    )
    + "\n"
)
