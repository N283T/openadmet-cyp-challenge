"""Count source-preserving external D6 measurements near known hit anchors."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    assignments = pl.read_parquet(inp / "parent-assignments.parquet")
    challenge = pl.read_parquet(inp / "challenge-compounds.parquet")
    anchor_ids = assignments.filter(pl.col("is_anchor")).select("compound_id").unique()
    anchors = anchor_ids.join(challenge, on="compound_id", validate="1:1")
    if anchors.height != 70:
        raise ValueError("Unexpected challenge anchor count")
    observations = pl.read_parquet(inp / "external-observations.parquet").filter(
        (pl.col("isoform") == cfg["isoform"])
        & pl.col("source").is_in(cfg["sources"])
        & pl.col("measurement_type").is_in(cfg["measurement_types"])
        & (
            pl.col("value").is_not_null()
            | pl.col("lower_bound").is_not_null()
            | pl.col("upper_bound").is_not_null()
        )
    ).select("compound_id", "connectivity", "challenge_overlap", "observation_id",
             "source", "assay_id", "measurement_type", "value", "relation",
             "lower_bound", "upper_bound")
    if observations.filter(pl.col("challenge_overlap")).height:
        raise ValueError("Challenge-overlap observations in external corpus")
    external = pl.read_parquet(inp / "external-compounds.parquet")
    measured = observations.select("compound_id").unique().join(
        external.select("compound_id", "standardized_smiles", "connectivity"),
        on="compound_id", validate="1:1"
    ).sort("compound_id")
    if measured.height != observations["compound_id"].n_unique():
        raise ValueError("Missing measured external compounds")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    anchor_fps = []
    for smiles in anchors["standardized_smiles"]:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError("Invalid anchor SMILES")
        anchor_fps.append(generator.GetFingerprint(molecule))
    nearest = []
    for row in measured.iter_rows(named=True):
        molecule = Chem.MolFromSmiles(row["standardized_smiles"])
        if molecule is None:
            raise ValueError(f"Invalid external SMILES {row['compound_id']}")
        fp = generator.GetFingerprint(molecule)
        similarities = DataStructs.BulkTanimotoSimilarity(fp, anchor_fps)
        best = int(np.argmax(similarities))
        nearest.append({
            "compound_id": row["compound_id"],
            "nearest_anchor_id": anchors["compound_id"][best],
            "max_tanimoto": float(similarities[best]),
        })
    proximity = pl.DataFrame(nearest)
    annotated = observations.join(proximity, on="compound_id", validate="m:1")
    if annotated.height != observations.height:
        raise ValueError("Incomplete observation proximity join")
    rows = []
    for (source, measurement_type), frame in annotated.group_by("source", "measurement_type"):
        for threshold in [0.0] + cfg["tanimoto_thresholds"]:
            part = frame.filter(pl.col("max_tanimoto") >= threshold)
            exact = part.filter(pl.col("relation") == "=")
            weak = exact.filter(pl.col("value") < cfg["weak_pic50_cutoff"])
            censored_below_4 = part.filter(
                pl.col("relation").is_in(["<", "<="])
                & (pl.col("upper_bound") <= cfg["weak_pic50_cutoff"])
            )
            censored_below_5 = part.filter(
                pl.col("relation").is_in(["<", "<="])
                & (pl.col("upper_bound") <= 5.0)
            )
            rows.append({
                "source": source, "measurement_type": measurement_type,
                "minimum_tanimoto": threshold,
                "observations": part.height,
                "compounds": part["compound_id"].n_unique(),
                "assays": part["assay_id"].n_unique(),
                "anchors_reached": part["nearest_anchor_id"].n_unique(),
                "exact_observations": exact.height,
                "exact_compounds": exact["compound_id"].n_unique(),
                "weak_exact_observations": weak.height,
                "weak_exact_compounds": weak["compound_id"].n_unique(),
                "censored_below_4_observations": censored_below_4.height,
                "censored_below_4_compounds": censored_below_4["compound_id"].n_unique(),
                "censored_below_5_observations": censored_below_5.height,
                "censored_below_5_compounds": censored_below_5["compound_id"].n_unique(),
                "exact_mean_value": exact["value"].mean(),
                "exact_median_value": exact["value"].median(),
            })
    pl.DataFrame(rows).sort("source", "measurement_type", "minimum_tanimoto").write_csv(
        out / "coverage.csv"
    )
    annotated.write_parquet(out / "observations-with-nearest-anchor.parquet")
    (out / "summary.json").write_text(json.dumps({
        "anchors": anchors.height,
        "measured_external_compounds": measured.height,
        "external_observations": observations.height,
        "fingerprint": cfg["fingerprint"],
        "scope": cfg["scope"],
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
