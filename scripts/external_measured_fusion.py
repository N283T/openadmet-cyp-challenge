"""External measured-value availability features, never challenge targets.

Keep the twenty external tree predictions and add measured-or-predicted values
plus measured-value availability flags. Match exact standardized identities only.
"""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    pl.Config.set_tbl_rows(25)

    def path(name):
        matches = [p for p in inp.iterdir() if p.name.partition("-")[2] == name]
        assert len(matches) == 1
        return matches[0]

    def read(name):
        return pl.read_parquet(path(name))

    compounds = read("compounds.parquet").sort("compound_id")
    test = read("test_mapping.parquet")
    predictions = read("external-scores.parquet").sort("compound_id")
    tree_metadata = json.loads(path("summary.json").read_text())
    assert predictions["compound_id"].equals(compounds["compound_id"])
    assert predictions["compound_id"].n_unique() == predictions.height
    assert len(tree_metadata["models"]) == 20
    pubchem = read("pubchem-curated-records.parquet").filter(
        pl.col("compound_id").is_in(compounds["compound_id"].implode())
    )
    tox21 = read("tox21-curated-records.parquet").filter(
        pl.col("compound_id").is_in(compounds["compound_id"].implode())
    )
    replicates = read("tox21-curated-replicates.parquet").filter(
        pl.col("compound_id").is_in(compounds["compound_id"].implode())
    )
    pubchem.write_parquet(out / "matched-pubchem-source-records.parquet")
    tox21.write_parquet(out / "matched-tox21-source-records.parquet")
    replicates.write_parquet(out / "matched-tox21-source-replicates.parquet")
    pub = pubchem.with_columns(
        (pl.lit("pubchem1851:") + pl.col("Panel ID") + ":" + pl.col("isoform") + ":AC50").alias(
            "assay_id"
        ),
        pl.when(pl.col("PUBCHEM_ACTIVITY_OUTCOME") != "Inactive")
        .then(-pl.col("Fit_LogAC50").cast(pl.Float64, strict=False))
        .alias("potency"),
        pl.when(pl.col("PUBCHEM_ACTIVITY_OUTCOME") == "Active")
        .then(1.0)
        .when(pl.col("PUBCHEM_ACTIVITY_OUTCOME") == "Inactive")
        .then(0.0)
        .alias("binary"),
    )
    tox = tox21.with_columns(
        (
            pl.lit("tox21:") + pl.col("aid").cast(pl.String) + ":" + pl.col("isoform") + ":AC50"
        ).alias("assay_id"),
        pl.when(pl.col("PUBCHEM_ACTIVITY_OUTCOME") == "Active")
        .then(1.0)
        .when(pl.col("PUBCHEM_ACTIVITY_OUTCOME") == "Inactive")
        .then(0.0)
        .alias("binary"),
    )
    active = tox21.filter(pl.col("PUBCHEM_ACTIVITY_OUTCOME") == "Active").select(
        "aid", "PUBCHEM_SID"
    )
    potency = (
        replicates.join(active, on=["aid", "PUBCHEM_SID"], validate="m:1")
        .filter(
            (pl.col("Phenotype") == "Inhibitor")
            & pl.col("pac50").is_finite()
            & pl.col("curve_class").is_in([-1.1, -1.2, -2.1, -2.2])
        )
        .with_columns(
            (
                pl.lit("tox21:") + pl.col("aid").cast(pl.String) + ":" + pl.col("isoform") + ":AC50"
            ).alias("assay_id"),
            pl.col("pac50").alias("potency"),
        )
    )
    observed, coverage = [], []
    test_ids = set(test["compound_id"])
    role = np.array(["test" if key in test_ids else "non_test" for key in compounds["compound_id"]])
    result = predictions
    for model in tree_metadata["models"]:
        name, assay, source, objective = (
            model["name"],
            model["assay_id"],
            model["source"],
            model["objective"],
        )
        if objective == "binary":
            raw = pub if source == "pubchem1851" else tox
            selected = (
                raw.filter(pl.col("assay_id") == assay)
                .group_by("compound_id")
                .agg(
                    pl.when(pl.col("binary").drop_nulls().n_unique() == 1)
                    .then(pl.col("binary").drop_nulls().first())
                    .otherwise(None)
                    .alias("measured"),
                    (pl.col("binary").drop_nulls().n_unique() > 1).alias("conflict"),
                    pl.col("PUBCHEM_SID").cast(pl.String).unique().alias("source_sids"),
                    pl.len().alias("source_records"),
                )
            )
        else:
            raw = pub if source == "pubchem1851" else potency
            selected = (
                raw.filter((pl.col("assay_id") == assay) & pl.col("potency").is_finite())
                .group_by("compound_id")
                .agg(
                    pl.col("potency").median().alias("measured"),
                    pl.lit(False).alias("conflict"),
                    pl.col("PUBCHEM_SID").cast(pl.String).unique().alias("source_sids"),
                    pl.len().alias("source_records"),
                )
            )
        selected = selected.with_columns(
            pl.lit(name).alias("score_name"),
            pl.lit(assay).alias("assay_id"),
            pl.lit(source).alias("source"),
            pl.lit(objective).alias("objective"),
        )
        observed.append(selected)
        joined = compounds.select("compound_id").join(
            selected.select("compound_id", "measured"),
            on="compound_id",
            how="left",
            validate="1:1",
            maintain_order="left",
        )
        measured = joined["measured"].to_numpy()
        available = np.isfinite(measured)
        original = predictions[name].to_numpy()
        result = result.with_columns(
            pl.Series(
                name.replace("score_", "score_fused_", 1), np.where(available, measured, original)
            ),
            pl.Series(name.replace("score_", "score_available_", 1), available.astype(np.int8)),
        )
        for subset in ["all", "non_test", "test"]:
            mask = np.ones(len(role), dtype=bool) if subset == "all" else role == subset
            coverage.append(
                dict(
                    score_name=name,
                    assay_id=assay,
                    source=source,
                    objective=objective,
                    subset=subset,
                    compounds=int(mask.sum()),
                    available=int((available & mask).sum()),
                    binary_conflicts=selected.filter("conflict").height
                    if subset == "all"
                    else None,
                )
            )
    assert result.width == 61 and result.height == 6896
    assert np.isfinite(result.drop("compound_id").to_numpy()).all()
    assert result.select(predictions.columns).equals(predictions)
    result.write_parquet(out / "external-measured-features.parquet")
    pl.concat(observed, how="vertical_relaxed").write_parquet(
        out / "observed-external-values.parquet"
    )
    pl.DataFrame(coverage).write_csv(out / "availability.csv")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                rows=result.height,
                original_predictions=20,
                fused_features=20,
                availability_flags=20,
                challenge_labels_used=False,
                matching="Exact standardized compound_id; no connectivity substitution",
                upstream_models="Challenge connectivity excluded in tree training; no retraining",
                policy="Same QC and aggregation for train/test; no fitted calibration",
                pubchem_policy="Finite -Fit_LogAC50 except Inactive; confident binary only",
                tox21_policy="Active Inhibitor quality-curve median; conflicting binary masked",
                coverage=coverage,
            ),
            indent=2,
        )
    )
    print(pl.DataFrame(coverage).filter(pl.col("subset") == "all"), flush=True)


if __name__ == "__main__":
    main()
