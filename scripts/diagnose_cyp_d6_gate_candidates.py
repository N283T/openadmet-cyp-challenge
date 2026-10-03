"""Diagnose D6 calls added or removed by the prepared gate candidates versus Sep 15.

For each prepared candidate the blind D6 calls are classified as retained, removed
or added relative to the Sep 15 payload. Each changed call is listed with its
legacy, pool_lgbm and combined ranks, the Morgan-gate p_active of run 7920b49a and
the similarity-weighted top-10 positive rate among labelled D6 training molecules.
Group means summarize retained, removed and added sets. Read-only: nothing is
fitted, registered, prepared or uploaded.
"""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")
GROUP_METRICS = [
    "legacy_rank",
    "pool_rank",
    "combined_rank",
    "legacy_position",
    "pool_position",
    "combined_position",
    "p_active_morgan",
    "top10_weighted_positive_rate",
    "nearest_positive_similarity",
]


def percentile_rank(values):
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values))
    ranks[order] = np.arange(len(values))
    return (ranks + 0.5) / len(values)


def positions(score):
    return np.argsort(np.argsort(-score, kind="stable"), kind="stable") + 1


def fingerprints(smiles, generator):
    fps = []
    for value in smiles:
        mol = Chem.MolFromSmiles(value)
        if mol is None:
            raise ValueError(f"Unparseable standardized SMILES: {value}")
        fps.append(generator.GetFingerprint(mol))
    return fps


def neighbour_table(inp, cfg, test_ids):
    compounds = pl.read_parquet(inp / "compounds.parquet")
    labels = (
        pl.read_parquet(inp / "measurements.parquet")
        .filter(pl.col("endpoint") == cfg["endpoint"])
        .select("compound_id", pl.col("value").cast(pl.Int8).alias("y"))
    )
    assert labels["compound_id"].n_unique() == labels.height
    assert set(labels["y"].unique().to_list()) == {0, 1}
    train = labels.join(compounds, on="compound_id", validate="1:1").sort("compound_id")
    assert train.height == labels.height
    test = test_ids.join(compounds, on="compound_id", how="left", validate="1:1")
    assert test["standardized_smiles"].null_count() == 0
    assert not set(test["compound_id"]) & set(train["compound_id"])
    fp = cfg["fingerprint"]
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=fp["radius"], fpSize=fp["fp_size"])
    train_fps = fingerprints(train["standardized_smiles"], generator)
    test_fps = fingerprints(test["standardized_smiles"], generator)
    y = train["y"].to_numpy()
    k = cfg["knn_k"]
    rows = []
    for compound_id, query in zip(test["compound_id"], test_fps, strict=True):
        sims = np.asarray(DataStructs.BulkTanimotoSimilarity(query, train_fps))
        top = np.argsort(-sims, kind="stable")[:k]
        weights = sims[top]
        rows.append(
            dict(
                compound_id=compound_id,
                top10_weighted_positive_rate=float(np.average(y[top], weights=weights))
                if weights.sum() > 0
                else float(y[top].mean()),
                top10_positive_rate=float(y[top].mean()),
                nearest_any_similarity=float(sims[top[0]]),
                nearest_positive_similarity=float(sims[y == 1].max()),
            )
        )
    meta = dict(
        labelled_molecules=train.height,
        labelled_positives=int(y.sum()),
        test_molecules=test.height,
    )
    return pl.DataFrame(rows), meta


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    endpoint = cfg["endpoint"]
    sep15 = pl.read_parquet(inp / "sep15_payload.parquet").select(
        pl.col("Molecule_Name").alias("molecule_name"), pl.col(endpoint).alias("sep15_payload")
    )
    gate_scores = pl.read_parquet(inp / "test_scores.parquet")
    gate = gate_scores.filter(pl.col("model") == cfg["gate_model"]).select(
        "compound_id", pl.col(cfg["gate_field"]).alias("gate_p_active")
    )
    pool = gate_scores.filter(pl.col("model") == cfg["pool_model"]).select(
        "compound_id", pl.col("score").alias("gate_run_pool_score")
    )
    assert gate.height == 750 and pool.height == 750

    neighbours, neighbour_meta = None, None
    tables, means, checks = [], [], {}
    for variant, spec in cfg["candidates"].items():
        scores = pl.read_parquet(inp / spec["scores_input"]).rename(
            {"Molecule_Name": "molecule_name"}
        )
        assert scores.height == 750 and scores["model"].unique().to_list() == [variant]
        if neighbours is None:
            neighbours, neighbour_meta = neighbour_table(inp, cfg, scores.select("compound_id"))
        prepared = pl.read_parquet(inp / spec["prepared_input"]).select(
            pl.col("Molecule_Name").alias("molecule_name"), pl.col(endpoint).alias("prepared")
        )
        frame = (
            scores.join(prepared, on="molecule_name", validate="1:1")
            .join(sep15, on="molecule_name", validate="1:1")
            .join(gate, on="compound_id", validate="1:1")
            .join(pool, on="compound_id", validate="1:1")
            .join(neighbours, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        assert frame.height == 750
        # The packaging-run scores must agree with the prepared file and the raw inputs.
        assert frame["call"].equals(frame["prepared"]), variant
        assert frame["sep15_call"].equals(frame["sep15_payload"]), variant
        assert frame["p_active_morgan"].equals(frame["gate_p_active"]), variant
        assert frame["pool_score"].equals(frame["gate_run_pool_score"]), variant
        score = frame["score"].to_numpy()
        frame = frame.with_columns(
            pl.Series("combined_rank", percentile_rank(score)),
            pl.Series("pool_position", positions(frame["pool_rank"].to_numpy())),
            pl.when(pl.col("call") & pl.col("sep15_call"))
            .then(pl.lit("retained"))
            .when(pl.col("sep15_call"))
            .then(pl.lit("removed"))
            .when(pl.col("call"))
            .then(pl.lit("added"))
            .otherwise(None)
            .alias("group"),
        ).rename({"position": "combined_position"})
        assert np.array_equal(frame["combined_position"].to_numpy(), positions(score))
        groups = frame.filter(pl.col("group").is_not_null()).with_columns(
            pl.lit(variant).alias("variant"),
            pl.lit(spec["prepared_id"]).alias("prepared_id"),
            pl.lit(endpoint).alias("endpoint"),
            pl.lit("blind").alias("split_id"),
        )
        counts = groups["group"].value_counts().sort("group")
        checks[variant] = dict(zip(counts["group"], counts["count"], strict=True))
        tables.append(groups)
        means.append(
            groups.group_by("group")
            .agg(pl.len().alias("n"), *[pl.col(c).mean() for c in GROUP_METRICS])
            .with_columns(pl.lit(variant).alias("variant"))
        )

    columns = [
        "variant",
        "prepared_id",
        "endpoint",
        "split_id",
        "group",
        "molecule_name",
        "compound_id",
        "legacy_rank",
        "legacy_position",
        "pool_rank",
        "pool_position",
        "combined_rank",
        "combined_position",
        "p_active_morgan",
        "top10_weighted_positive_rate",
        "top10_positive_rate",
        "nearest_any_similarity",
        "nearest_positive_similarity",
    ]
    order = {"retained": 0, "removed": 1, "added": 2}
    calls = (
        pl.concat(tables)
        .select(columns)
        .with_columns(pl.col("group").replace_strict(order).alias("_order"))
        .sort("variant", "_order", "combined_position", descending=[True, False, False])
        .drop("_order")
    )
    calls.write_parquet(out / "call_diagnostics.parquet")
    calls.filter(pl.col("group") != "retained").write_csv(out / "changed_calls.csv")
    group_means = (
        pl.concat(means)
        .with_columns(pl.col("group").replace_strict(order).alias("_order"))
        .sort("variant", "_order", descending=[True, False])
        .drop("_order")
        .select("variant", "group", "n", *GROUP_METRICS)
    )
    group_means.write_csv(out / "group_means.csv")
    neighbours.write_parquet(out / "blind_neighbour_rates.parquet")
    summary = dict(
        group_counts=checks,
        neighbour_reference=neighbour_meta,
        group_means=group_means.to_dicts(),
        uploaded=False,
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
