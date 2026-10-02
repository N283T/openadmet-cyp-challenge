"""Check the refit composition against archived full-pool atom predictions."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from validate_cyp_d6_legacy_ensemble import compose_legacy


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    sources = {
        name: pl.read_parquet(inp / f"{name}.parquet")
        for name in ["lgbm", "pfn_pca", "pfn_raw", "graph", "fixed", "control"]
    }
    members = [
        ("lgbm", m)
        for m in [
            "morgan_plain",
            "mordred_plain",
            "fp_mordred_plain",
            "external_class_plain",
            "external_mixed_plain",
        ]
    ]
    members += [
        ("pfn_pca", "external_mixed_pca256"),
        ("pfn_raw", "external_mixed_raw600"),
        ("graph", "tdi_only"),
        ("fixed", "logistic_c0.05"),
        ("fixed", "rbf_svm_c1"),
        ("fixed", "complement_nb_a1"),
    ]
    parts = [
        sources[source].filter(pl.col("model") == model).sort("compound_id")
        for source, model in members
    ]
    reference = (
        sources["control"]
        .filter(pl.col("model") == "sep15_control_raw_oof_analog")
        .sort("compound_id")
    )
    keys = ["compound_id", "endpoint", "split_id", "fold", "group"]
    for part in parts:
        assert part.height == 6146 and part.select(keys).equals(reference.select(keys))
    raw = np.column_stack([part["score"].to_numpy() for part in parts])
    actual = np.full(reference.height, np.nan)
    for fold in range(5):
        mask = reference["fold"].to_numpy() == fold
        actual[mask] = compose_legacy(raw[mask])
    expected = reference["score"].to_numpy()
    error = float(np.max(np.abs(actual - expected)))
    assert error < 1e-12
    reference.select(keys).with_columns(
        pl.Series("expected", expected), pl.Series("actual", actual)
    ).write_parquet(out / "composition-audit.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                rows=6146,
                atoms=11,
                max_absolute_error=error,
                exact_composition_verified=True,
                fitting_or_decision_claim=False,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
