"""Blind D6 calls from the legacy Sep 15 score and a MONROE-pretrained pool LightGBM.

The deployed L4 score is the equal percentile-rank mean of the legacy Sep 15 blind score and the
`pool_lgbm` blind score. Here `pool_lgbm` is replaced by a balanced LightGBM on frozen MONROE
vectors (after external CYP pretraining), fitted on every labelled pool row (official labels plus
primary-screen pseudo-negatives). The L4 construction is first reproduced against the last sent
payload. A4 is copied unchanged. Nothing is uploaded.
"""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from lightgbm import LGBMClassifier
from package_cyp_d6_gate_candidate import (
    IDENTIFIERS,
    blind_column,
    payload,
    percentile_rank,
    top_calls,
)

from cyp_challenge.tracks import validate


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    d6, a4 = cfg["d6_endpoint"], cfg["a4_endpoint"]
    reference = pl.read_csv(inp / "reference.csv").select(IDENTIFIERS).sort("Molecule_Name")
    names = (
        pl.read_parquet(inp / "test_mapping.parquet")
        .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
        .unique()
    )
    last = payload(inp / "last_sent.parquet", reference)
    assert int(last[d6].sum()) == cfg["calls"] and int(last[a4].sum()) == cfg["a4_positive_calls"]

    legacy = pl.read_parquet(inp / "sep15_blind_scores.parquet").filter(
        pl.col("model") == cfg["legacy_model"]
    )
    table = legacy.select("compound_id", pl.col("score").alias("legacy_score")).join(
        names, on="compound_id", how="left", validate="1:1"
    )
    table = table.join(
        blind_column(pl.read_parquet(inp / "test_scores.parquet"), cfg["pool_column"]),
        on="compound_id",
        validate="1:1",
    ).sort("compound_id")
    assert table.height == 750
    legacy_rank = percentile_rank(table["legacy_score"].to_numpy())
    pool_rank = percentile_rank(table[cfg["pool_column"]].to_numpy())
    last_call = (
        table.select("Molecule_Name")
        .join(last.select("Molecule_Name", pl.col(d6)), on="Molecule_Name", maintain_order="left")[
            d6
        ]
        .to_numpy()
    )
    l4 = 0.5 * (legacy_rank + pool_rank)
    assert np.array_equal(top_calls(l4, cfg["calls"]), last_call), "L4 not reproduced"

    labels = (
        pl.read_parquet(inp / "pool-labels.parquet")
        .filter((pl.col("split_id") == cfg["label_split"]) & (pl.col("model") == "pool_lgbm"))
        .select("compound_id", "y_testlike")
        .filter(pl.col("y_testlike").is_not_null())
    )
    embeddings = pl.read_parquet(inp / "embeddings.parquet")
    columns = [c for c in embeddings.columns if c.startswith("emb_")]
    train = labels.join(embeddings, on="compound_id", how="inner", validate="1:1").sort(
        "compound_id"
    )
    assert train.height == labels.height
    x_blind = table.select("compound_id").join(embeddings, on="compound_id", validate="1:1")
    scores = []
    for s in range(cfg["lgbm_seeds"]):
        model = LGBMClassifier(
            class_weight="balanced",
            **cfg["lgbm"],
            verbosity=-1,
            n_jobs=cfg["threads"],
            random_state=cfg["seed"] + s,
            deterministic=True,
            force_col_wise=True,
        )
        model.fit(
            train.select(columns).to_numpy().astype(np.float32),
            train["y_testlike"].to_numpy().astype(int),
        )
        scores.append(
            np.asarray(model.predict_proba(x_blind.select(columns).to_numpy().astype(np.float32)))[
                :, 1
            ]
        )
    monroe = np.mean(scores, axis=0)
    swap = 0.5 * (legacy_rank + percentile_rank(monroe))
    count = cfg["calls"]
    calls = top_calls(swap, count)
    order = np.argsort(-swap, kind="stable")
    gap = float(swap[order[count - 1]] - swap[order[count]])
    assert gap > 0, "Tied call boundary"

    candidate = (
        last.drop(d6)
        .join(
            table.select("Molecule_Name").with_columns(pl.Series(d6, calls)),
            on="Molecule_Name",
            validate="1:1",
            maintain_order="left",
        )
        .select(*IDENTIFIERS, d6, a4)
    )
    assert candidate[a4].equals(last[a4]) and int(candidate[d6].sum()) == count
    candidate.write_parquet(out / "predictions.parquet")
    candidate.write_csv(out / "predictions.csv")
    validation = validate(out / "predictions.parquet", "classification", inp / "reference.csv")
    summary = {
        "calls": count,
        "kept_from_last_sent": int((calls & last_call).sum()),
        "added": int((calls & ~last_call).sum()),
        "removed": int((~calls & last_call).sum()),
        "boundary_gap": gap,
        "train_rows": train.height,
        "train_positives": int(train["y_testlike"].sum()),
        "spearman_legacy_vs_monroe": float(
            pl.DataFrame({"a": legacy_rank, "b": percentile_rank(monroe)})
            .select(pl.corr("a", "b", method="spearman"))
            .item()
        ),
        "spearman_pool_vs_monroe": float(
            pl.DataFrame({"a": pool_rank, "b": percentile_rank(monroe)})
            .select(pl.corr("a", "b", method="spearman"))
            .item()
        ),
        "a4_unchanged": True,
        "validation": validation,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps(summary, default=str))


if __name__ == "__main__":
    main()
