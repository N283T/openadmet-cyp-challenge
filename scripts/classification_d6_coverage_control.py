"""Match test neighbor pool to heldout UMAP pools; distinguish recall and FPR by assay route."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


comp, labels, assign, test, route, oof = [pl.read_parquet(path(i)) for i in range(6)]
assign = assign.filter(pl.col("strategy") == "umap")
comp = comp.sort("compound_id")
ids = comp["compound_id"].to_numpy()
lookup = {c: i for i, c in enumerate(ids)}
fg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps = [fg.GetFingerprint(Chem.MolFromSmiles(s)) for s in comp["standardized_smiles"]]
testix = [lookup[c] for c in test["compound_id"]]
coverage = []
pairs = []
for ep in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]:
    labelled = labels.filter(pl.col("endpoint") == ep).join(
        assign, on="compound_id", validate="1:1"
    )
    for fold in range(5):
        pool = [lookup[c] for c in labelled.filter(pl.col("fold") != fold)["compound_id"]]
        poolfp = [fps[i] for i in pool]
        for role, query in [
            ("test", testix),
            (
                "heldout",
                [lookup[c] for c in labelled.filter(pl.col("fold") == fold)["compound_id"]],
            ),
        ]:
            sim = np.array([max(DataStructs.BulkTanimotoSimilarity(fps[i], poolfp)) for i in query])
            coverage.append(
                dict(
                    endpoint=ep,
                    fold=fold,
                    role=role,
                    n=len(query),
                    pool_size=len(pool),
                    median_similarity=float(np.median(sim)),
                    mean_similarity=float(sim.mean()),
                    fraction_ge05=float((sim >= 0.5).mean()),
                )
            )
            pairs.append(
                pl.DataFrame(
                    dict(compound_id=ids[query], endpoint=ep, fold=fold, role=role, similarity=sim)
                )
            )
pl.DataFrame(coverage).write_csv(out / "matched-pool-coverage.csv")
pl.concat(pairs).write_parquet(out / "matched-pool-neighbors.parquet")
errors = []
for ep in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]:
    method = "all_plain_pfn_graph_mean" if "2D6" in ep else "lgbm_balanced_mean"
    h = oof.filter((pl.col("endpoint") == ep) & (pl.col("method") == method)).join(
        route.filter(pl.col("endpoint") == ep).select(
            "compound_id", "availability", "shift_regime", "bounds_regime"
        ),
        on="compound_id",
        validate="1:1",
    )
    for dimension in ["availability", "shift_regime", "bounds_regime"]:
        for (key,), t in h.group_by(dimension):
            y = t["y_true"].to_numpy().astype(bool)
            p = t["prediction"].to_numpy() >= t["threshold"].to_numpy()
            errors.append(
                dict(
                    endpoint=ep,
                    dimension=dimension,
                    stratum=key,
                    n=len(t),
                    positive=int(y.sum()),
                    negative=int((~y).sum()),
                    TP=int((p & y).sum()),
                    FN=int((~p & y).sum()),
                    FP=int((p & ~y).sum()),
                    TN=int((~p & ~y).sum()),
                    recall=float(p[y].mean()) if y.any() else None,
                    FPR=float(p[~y].mean()) if (~y).any() else None,
                )
            )
pl.DataFrame(errors).write_csv(out / "route-confusion-counts.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            test_pool_control="Every fold uses exactly same labelled training pool for heldout "
            "and test; five test measurements repeated, not independent",
            canonical_missingness="Missing canonical measurement does not establish that "
            "assay was not performed",
        ),
        indent=2,
    )
)
