"""Change only the submitted A4 equal-seven decision count to 275."""

import json
import os
from pathlib import Path

import polars as pl

from cyp_challenge.tracks import ENDPOINTS, validate


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    d6, a4 = ENDPOINTS["classification"]
    current = pl.read_parquet(inp / cfg["current_candidate"]).sort("Molecule_Name")
    scores = pl.read_parquet(inp / cfg["submitted_scores"])
    aliases = (
        pl.read_parquet(inp / cfg["aliases"])
        .filter(pl.col("role") == "test")
        .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
        .unique()
    )
    if current.height != 750 or scores.height != 750 or aliases.height != 750:
        raise ValueError("Expected 750 distinct blind molecules")
    if scores["compound_id"].n_unique() != 750 or aliases["Molecule_Name"].n_unique() != 750:
        raise ValueError("Duplicate blind identity")
    ranked = scores.sort(["a4_prediction", "compound_id"], descending=[True, False])
    top_reference = ranked.head(cfg["reference_a4_calls"]).select("compound_id")
    original = scores.select("compound_id", "a4_call").join(
        top_reference.with_columns(pl.lit(True).alias("replayed")),
        on="compound_id", how="left", validate="1:1"
    ).with_columns(pl.col("replayed").fill_null(False))
    if not original["a4_call"].equals(original["replayed"]):
        raise ValueError("Submitted A4 calls are not the deterministic top-249 scores")
    selected = ranked.head(cfg["a4_positive_calls"]).select("compound_id")
    calls = aliases.join(
        selected.with_columns(pl.lit(True).alias("new_a4")),
        on="compound_id", how="left", validate="1:1"
    ).with_columns(pl.col("new_a4").fill_null(False)).select("Molecule_Name", "new_a4")
    candidate = current.join(calls, on="Molecule_Name", validate="1:1").drop(a4).rename(
        {"new_a4": a4}
    ).select("SMILES", "Molecule_Name", d6, a4).sort("Molecule_Name")
    if not candidate[d6].equals(current[d6]):
        raise ValueError("D6 changed")
    if candidate[a4].sum() != cfg["a4_positive_calls"]:
        raise ValueError("Wrong A4 call count")
    if ((current[a4]) & (~candidate[a4])).sum() != 0:
        raise ValueError("A4 rate candidate removed a submitted positive")
    curve = pl.read_csv(inp / cfg["rate_curve"])
    target = curve.filter(
        (pl.col("prior") == 0.2916) & (pl.col("count_per_750") == cfg["a4_positive_calls"])
    )
    if target.height != 1 or target["delta_vs_reference"][0] <= 0:
        raise ValueError("Selected rate lacks prespecified OOF support")
    candidate.write_parquet(out / "predictions.parquet")
    candidate.write_csv(out / "predictions.csv")
    summary = {
        "selection": cfg["selection"],
        "positive_calls": {d6: int(candidate[d6].sum()), a4: int(candidate[a4].sum())},
        "new_a4_calls": int((~current[a4] & candidate[a4]).sum()),
        "removed_a4_calls": int((current[a4] & ~candidate[a4]).sum()),
        "d6_bitwise_unchanged": True,
        "submitted_top_249_reproduced": True,
        "oof_delta_mcc_at_working_prior": float(target["delta_vs_reference"][0]),
        "validation": validate(
            out / "predictions.parquet", "classification", inp / cfg["reference"]
        ),
        "uploaded": False,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
