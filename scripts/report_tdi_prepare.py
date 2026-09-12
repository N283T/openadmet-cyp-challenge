"""Extract paired continuous TDI auxiliaries with explicit source identity."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])

    def read(name):
        paths = list(inp.glob("*-" + name))
        assert len(paths) == 1
        return paths[0]

    raw = pl.read_csv(read("cyp-challenge-TRAIN_TDI.csv"))
    aliases = (
        pl.read_parquet(read("aliases.parquet"))
        .select(
            pl.col("molecule_name").alias("Molecule_Name"),
            pl.col("raw_smiles").alias("SMILES"),
            "compound_id",
        )
        .unique()
    )
    raw = raw.join(aliases, on=["Molecule_Name", "SMILES"], how="left", validate="m:1")
    assert raw["compound_id"].null_count() == 0
    labels = pl.read_parquet(read("measurements.parquet"))
    assignments = pl.read_parquet(read("parent-assignments.parquet"))
    arms, deltas, summary = [], [], []
    for target in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]:
        direct, tdi = target + "_pIC50_direct_inhibition", target + "_pIC50_TDI_condition"
        part = raw.select("compound_id", direct, tdi).unique()
        assert part["compound_id"].n_unique() == part.height
        check = part.filter(pl.col(direct).is_finite()).join(
            labels.filter(pl.col("endpoint") == direct).select("compound_id", "value"),
            on="compound_id",
            how="left",
            validate="1:1",
        )
        assert check["value"].null_count() == 0
        assert np.allclose(check[direct], check["value"], atol=1e-10, rtol=0)
        arm = part.filter(pl.col(tdi).is_finite()).select(
            "compound_id", pl.lit(tdi).alias("head"), pl.col(tdi).alias("value")
        )
        paired = part.filter(pl.col(tdi).is_finite() & pl.col(direct).is_finite())
        delta = paired.select(
            "compound_id",
            pl.lit(target + "_pIC50_TDI_minus_direct").alias("head"),
            (pl.col(tdi) - pl.col(direct)).alias("value"),
        )
        arms.append(arm)
        deltas.append(delta)
        summary.append(
            dict(
                target=target,
                tdi_labels=arm.height,
                paired_labels=delta.height,
                tdi_only=arm.height - delta.height,
                direct_tdi_pearson=float(np.corrcoef(paired[direct], paired[tdi])[0, 1]),
            )
        )
    arm_table, delta_table = pl.concat(arms), pl.concat(deltas)
    for table in [arm_table, delta_table]:
        assert set(table["compound_id"]) <= set(assignments["compound_id"])
        assert table.select("compound_id", "head").n_unique() == table.height
    arm_table.write_parquet(out / "tdi-arm-targets.parquet")
    pl.concat([arm_table, delta_table]).write_parquet(out / "tdi-arm-delta-targets.parquet")
    raw.write_parquet(out / "mapped-tdi-source.parquet")
    result = dict(
        targets=summary,
        arm_labels=arm_table.height,
        delta_labels=delta_table.height,
        interpretation=(
            "TDI condition is +NADPH preincubation; direct is the supplied -NADPH comparator. "
            "Delta = TDI pIC50 minus direct pIC50. "
            "No inferred paired values or classification labels."
        ),
    )
    (out / "summary.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
