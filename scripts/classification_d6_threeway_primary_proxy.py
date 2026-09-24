"""Audit an independently trained primary-screen proxy on three D6 decisions."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
test = pl.read_parquet(inp / cfg["threeway_test"])
oof = pl.read_parquet(inp / cfg["threeway_oof"])
primary = pl.read_parquet(inp / cfg["primary_predictions"])

# For each labelled validation molecule use the auxiliary fit that excluded
# its outer fold. The blind set is scored by the mean of the five such fits.
fold_safe = primary.filter(pl.col("fold") == pl.col("outer_fold")).select(
    "compound_id",
    pl.col("fold").alias("primary_fold"),
    "is_anchor",
    pl.col(cfg["primary_column"]).alias("predicted_primary_log2fc"),
)
oof = oof.join(
    fold_safe,
    on="compound_id",
    validate="1:1",
)
if oof.height != 1497 or oof["predicted_primary_log2fc"].null_count():
    raise ValueError("Incomplete fold-excluded primary predictions")
honest_oof = oof.filter(~pl.col("is_anchor"))
if honest_oof.height < 1000:
    raise ValueError("Too few nonanchor OOF molecules")
blind = (
    primary.filter(pl.col("fold").is_null())
    .group_by("compound_id")
    .agg(
        pl.col("outer_fold").n_unique().alias("n_outer"),
        pl.col(cfg["primary_column"]).mean().alias("predicted_primary_log2fc"),
    )
)
if blind.height != 750 or not (blind["n_outer"] == 5).all():
    raise ValueError("Incomplete five-view blind predictions")
test = test.join(
    blind.select("compound_id", "predicted_primary_log2fc"), on="compound_id", validate="1:1"
)
if test.height != 750 or test["predicted_primary_log2fc"].null_count():
    raise ValueError("Incomplete blind join")

oof.write_parquet(out / "oof-primary.parquet")
test.write_parquet(out / "test-primary.parquet")

summary = {
    "nonanchor_oof_n": honest_oof.height,
    "excluded_anchor_oof_n": oof.height - honest_oof.height,
    "fold_excluded_primary_auc_for_D6_class": float(
        roc_auc_score(honest_oof["control_label"], -honest_oof["predicted_primary_log2fc"])
    ),
    "oof_positive_mean_log2fc": float(
        honest_oof.filter(pl.col("control_label") == 1)["predicted_primary_log2fc"].mean()
    ),
    "oof_negative_mean_log2fc": float(
        honest_oof.filter(pl.col("control_label") == 0)["predicted_primary_log2fc"].mean()
    ),
}

rng = np.random.default_rng(cfg["permutation_seed"])
contrasts = []
for population, table in [("blind", test), ("oof_nonanchor", honest_oof)]:
    for candidate in ["support", "ap"]:
        removed = table.filter(pl.col("control_call") & ~pl.col(f"{candidate}_call"))
        added = table.filter(~pl.col("control_call") & pl.col(f"{candidate}_call"))
        a = removed["predicted_primary_log2fc"].to_numpy()
        b = added["predicted_primary_log2fc"].to_numpy()
        observed = float(a.mean() - b.mean())
        pooled = np.concatenate([a, b])
        extreme = 0
        for _ in range(cfg["permutation_replicates"]):
            shuffle = rng.permutation(pooled)
            delta = shuffle[: len(a)].mean() - shuffle[len(a) :].mean()
            extreme += abs(delta) >= abs(observed)
        contrasts.append(
            {
                "population": population,
                "candidate": candidate,
                "removed_n": len(a),
                "added_n": len(b),
                "control_only_mean_log2fc": float(a.mean()),
                "candidate_only_mean_log2fc": float(b.mean()),
                "control_minus_candidate_log2fc": observed,
                "permutation_p_two_sided": (extreme + 1) / (cfg["permutation_replicates"] + 1),
            }
        )
pl.DataFrame(contrasts).write_csv(out / "primary-swap-contrasts.csv")
summary["interpretation"] = (
    "Lower predicted log2fc indicates greater predicted primary inhibition; "
    "this proxy is not a hidden class label."
)
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps({"summary": summary, "contrasts": contrasts}, indent=2), flush=True)
