"""Fixed seed/scope means and independent specialist seed checks."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    spec = importlib.util.spec_from_file_location(
        "official", next(inp.glob("*-custom_scoring_functions.py"))
    )
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    endpoint = "CYP2D6_pIC50_direct_inhibition"
    keys = ["compound_id", "endpoint", "fold", "split_id"]
    global_table = pl.read_parquet(inp / cfg["global_file"]).filter(pl.col("endpoint") == endpoint)
    a = pl.read_parquet(inp / cfg["assignments"])
    check = global_table.join(
        a.select("compound_id", pl.col("fold").alias("assigned_fold"), "is_anchor"),
        on="compound_id",
        validate="m:1",
    )
    assert check.height == global_table.height
    assert check.filter(
        (pl.col("fold") != pl.col("assigned_fold")) | pl.col("is_anchor")
    ).is_empty()
    global_table = global_table.with_columns(pl.lit(cfg["split_run"]).alias("split_id"))
    reference = global_table.filter(pl.col("method") == "caruana_weight_crossfit").sort(keys)
    tables = {}
    for name, filename in cfg["predictions"].items():
        table = pl.read_parquet(inp / filename).filter(pl.col("endpoint") == endpoint).sort(keys)
        assert table.height == 1457 and table.select(keys).equals(reference.select(keys))
        for col in ["y_true", "conf_low", "conf_high"]:
            assert np.array_equal(table[col].to_numpy(), reference[col].to_numpy())
        assert np.isfinite(table["prediction"].to_numpy()).all()
        tables[name] = table.with_columns(pl.lit(name).alias("method"))
    audit = []
    for seed, files in cfg["memberships"].items():
        for fold, filename in enumerate(files):
            m = pl.read_parquet(inp / filename).sort("compound_id")
            original = pl.read_parquet(inp / cfg["memberships"]["42"][fold]).sort("compound_id")
            assert m.equals(original), (seed, fold)
            outer = (pl.col("fold") == fold) & ~pl.col("is_anchor")
            assert m.filter(outer & (pl.col("role") != "outer_holdout")).is_empty()
            assert not set(m.filter(pl.col("role") == "training")["group"]) & set(
                m.filter(pl.col("role") == "inner_validation")["group"]
            )
            hist = pl.read_ndjson(inp / cfg["histories"][seed][fold])
            expected = hist["validation_pic50_head_st_rae"].struct.field(endpoint).to_numpy()
            assert np.array_equal(hist["selection_score"].to_numpy(), expected)
            audit.append(
                dict(
                    seed=int(seed),
                    fold=fold,
                    memberships_identical_to_seed42=True,
                    selection_2d6=True,
                )
            )
    for name, members in cfg["means"].items():
        values = np.mean([tables[member]["prediction"].to_numpy() for member in members], axis=0)
        table = reference.with_columns(
            pl.Series("prediction", values), pl.lit(name).alias("method")
        )
        table.write_parquet(out / (name + "-oof.parquet"))
        tables[name] = table
    metrics, rows = [], []
    for name, table in tables.items():
        for baseline_name in ["caruana_weight_crossfit", "caruana_full_fit"]:
            b = global_table.filter(pl.col("method") == baseline_name).sort(keys)
            assert table.select(keys).equals(b.select(keys))
            for weight in [0.0, 0.25, 0.5, 1.0]:
                prediction = (1 - weight) * b["prediction"].to_numpy() + weight * table[
                    "prediction"
                ].to_numpy()
                part = table.with_columns(pl.Series("prediction", prediction))
                rows.append(
                    part.select(
                        *keys, "y_true", "conf_low", "conf_high", "prediction", "method"
                    ).with_columns(
                        pl.lit(baseline_name).alias("baseline"), pl.lit(weight).alias("weight")
                    )
                )
                for fold in [-1, 0, 1, 2, 3, 4]:
                    f = part if fold == -1 else part.filter(pl.col("fold") == fold)
                    score = scorer.rae_soft_threshold_absolute_error(
                        f["y_true"].to_numpy(),
                        f["prediction"].to_numpy(),
                        y_true_lower=f["conf_low"].to_numpy(),
                        y_true_upper=f["conf_high"].to_numpy(),
                    )
                    metrics.append(
                        dict(
                            method=name,
                            baseline=baseline_name,
                            weight=weight,
                            fold=fold,
                            n=f.height,
                            ST_RAE=float(score),
                        )
                    )
    pl.DataFrame(metrics).write_csv(out / "metrics.csv")
    pl.concat(rows).write_parquet(out / "fixed-blends.parquet")
    (out / "audit.json").write_text(
        json.dumps(dict(checks=audit, fixed_means=cfg["means"], oof_rows=1457), indent=2)
    )
    print(json.dumps(dict(models=len(tables), audited_seed_folds=len(audit))), flush=True)


if __name__ == "__main__":
    main()
