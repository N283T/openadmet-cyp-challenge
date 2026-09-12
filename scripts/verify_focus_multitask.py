"""Verify persisted focus training memberships, histories, and OOF keys."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def main():
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def read_path(index):
        matches = list(inp.glob(f"{index:03d}-*"))
        assert len(matches) == 1
        return matches[0]

    summary = json.loads(read_path(0).read_text())
    config = summary["config"]
    endpoints = summary["heads"][:4]
    assert config["pic50_endpoint_weights"] == {
        "CYP1A2": 2.0,
        "CYP2C9": 0.5,
        "CYP2D6": 2.0,
        "CYP3A4": 0.5,
    }
    assert config["selection_targets"] == ["CYP1A2", "CYP2D6"]
    assert config["log2fc_weight"] == cfg["logfc_weight"]
    oof = pl.read_parquet(read_path(1))
    assert oof.height == 6384
    assert oof.select("compound_id", "endpoint").unique().height == 6384
    assert oof["prediction"].is_finite().all()
    records = []
    for fold in range(5):
        membership = pl.read_parquet(read_path(2 + fold * 2))
        history = [json.loads(line) for line in read_path(3 + fold * 2).read_text().splitlines()]
        heldout = membership.filter((pl.col("fold") == fold) & ~pl.col("is_anchor"))
        assert heldout["role"].unique().to_list() == ["outer_holdout"]
        fitted = membership.filter(pl.col("role").is_in(["training", "inner_validation"]))
        assert fitted.join(heldout, on="compound_id", how="inner").height == 0
        training_groups = set(membership.filter(pl.col("role") == "training")["group"])
        validation_groups = set(membership.filter(pl.col("role") == "inner_validation")["group"])
        assert not training_groups.intersection(validation_groups)
        predictions = oof.filter(pl.col("fold") == fold)
        assert (
            predictions.join(heldout.select("compound_id"), on="compound_id", how="anti").height
            == 0
        )
        assert set(predictions["endpoint"]) == set(endpoints)
        for row in history:
            scores = row["validation_pic50_head_st_rae"]
            expected = np.mean(
                [
                    value
                    for head, value in scores.items()
                    if head.split("_")[0] in config["selection_targets"]
                ]
            )
            assert abs(expected - row["selection_score"]) < 1e-12
            assert np.isfinite(row["train_pic50_interval"])
        best_epoch = min(history, key=lambda row: row["selection_score"])["epoch"]
        assert best_epoch == summary["folds"][fold]["best_epoch"]
        records.append(
            dict(
                fold=fold,
                oof_rows=predictions.height,
                best_epoch=best_epoch,
                epochs=len(history),
                outer_compounds=heldout.height,
            )
        )
    (out / "verification.json").write_text(
        json.dumps(dict(passed=True, rows=oof.height, config=config, folds=records), indent=2)
    )


if __name__ == "__main__":
    main()
