"""Independently audit canonical exports against native keyed all-head arrays."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    def read(name):
        return pl.read_parquet(inp / cfg[name])
    labels, parent, assignment = read("labels"), read("parent"), read("assignment")
    oof, test, fold_test = read("oof"), read("test"), read("fold_test")
    targets = ["CYP3A4", "CYP2C9", "CYP2D6", "CYP1A2"] if cfg["multi"] else ["CYP2D6"]
    endpoints = [target + "_pIC50_direct_inhibition" for target in targets]
    assert set(oof["endpoint"]) == set(test["endpoint"]) == set(endpoints)
    source_cv = cfg["source_cv"]
    assert (
        oof.height == (6525 if cfg["multi"] else 1493)
        if source_cv
        else oof.height == (6384 if cfg["multi"] else 1457)
    )
    assert oof.filter(pl.col("endpoint") == "CYP2D6_pIC50_direct_inhibition").height == (
        1493 if source_cv else 1457
    )
    assert oof.unique(["compound_id", "endpoint"]).height == oof.height
    assert test.height == 750 * len(targets)
    assert test.unique(["compound_id", "endpoint"]).height == test.height
    assert fold_test.height == 5 * test.height
    assert fold_test.unique(["compound_id", "endpoint", "fold"]).height == fold_test.height
    expected_test = set(read("mapping")["compound_id"])
    assert len(expected_test) == 750
    for endpoint in endpoints:
        assert set(test.filter(pl.col("endpoint") == endpoint)["compound_id"]) == expected_test
    expected = labels.filter(pl.col("endpoint").is_in(endpoints))
    if not source_cv:
        expected = expected.join(parent, on="compound_id", validate="m:1").filter(
            ~pl.col("is_anchor")
        )
        keys = ["compound_id", "endpoint", "fold"]
    else:
        expected = expected.join(
            assignment.select("compound_id", pl.col("effective_fold").alias("fold")),
            on="compound_id",
            validate="m:1",
        )
        keys = ["compound_id", "endpoint", "fold"]
    assert set(expected.select(keys).iter_rows()) == set(oof.select(keys).iter_rows())
    checked = oof.join(
        expected.select(
            *keys, "value", pl.col("conf_low").alias("low"), pl.col("conf_high").alias("high")
        ),
        on=keys,
        validate="1:1",
    )
    for actual, reference in [("y_true", "value"), ("conf_low", "low"), ("conf_high", "high")]:
        np.testing.assert_allclose(checked[actual], checked[reference], rtol=0, atol=0)
    deltas = []
    for fold in range(5):
        raw = read(f"val{fold}").join(
            assignment.select("row_id", "compound_id", "effective_fold"),
            on="row_id",
            validate="1:1",
        )
        assert (raw["effective_fold"] == fold).all()
        membership = read(f"member{fold}")
        assert set(raw["row_id"]) == set(membership.filter(pl.col("validation"))["row_id"])
        assert not membership.filter(pl.col("validation") & pl.col("training")).height
        raw_test = read(f"test{fold}")
        assert set(raw_test["compound_id"]) == expected_test
        for target, endpoint in zip(targets, endpoints, strict=True):
            head = target.lower() + "_pic50_direct_inhibition"
            actual = oof.filter((pl.col("fold") == fold) & (pl.col("endpoint") == endpoint))
            pair = actual.join(
                raw.select("compound_id", pl.col(head).alias("native")),
                on="compound_id",
                validate="1:1",
            )
            assert pair.height == actual.height and pair.height > 0
            np.testing.assert_array_equal(pair["prediction"].to_numpy(), pair["native"].to_numpy())
            other = fold_test.filter(
                (pl.col("fold") == fold) & (pl.col("endpoint") == endpoint)
            ).join(
                raw_test.select("compound_id", pl.col(head).alias("native")),
                on="compound_id",
                validate="1:1",
            )
            assert other.height == 750
            np.testing.assert_array_equal(
                other["prediction"].to_numpy(), other["native"].to_numpy()
            )
            deltas.append(
                dict(
                    fold=fold,
                    endpoint=endpoint,
                    oof_rows=pair.height,
                    test_rows=other.height,
                    max_absolute_difference=0.0,
                )
            )
    rebuilt = fold_test.group_by("compound_id", "endpoint").agg(
        pl.col("prediction").mean().alias("rebuilt"), pl.len().alias("count")
    )
    pair = test.join(rebuilt, on=["compound_id", "endpoint"], validate="1:1")
    assert (pair["count"] == 5).all()
    np.testing.assert_allclose(pair["prediction"], pair["rebuilt"], rtol=0, atol=1e-12)
    assert oof["prediction"].is_finite().all() and test["prediction"].is_finite().all()
    pl.DataFrame(deltas).write_csv(out / "mapping-checks.csv")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                training_run=cfg["training_run"],
                export_run=cfg["export_run"],
                oof_rows=oof.height,
                test_rows=test.height,
                raw_arrays_exact=True,
                canonical_membership_exact=True,
                labels_bounds_exact=True,
                test_mean_max_error=float((pair["prediction"] - pair["rebuilt"]).abs().max()),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
