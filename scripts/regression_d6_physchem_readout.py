"""Matched D6 direct-pIC50 readout screen with explicit molecular descriptors."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from catboost import CatBoostRegressor
from rdkit import Chem
from rdkit.Chem import Descriptors, Fragments, rdMolDescriptors
from scipy.stats import spearmanr

AMINE_SMARTS = {
    "phys_amine_primary": "[NX3;H2;!$(NC=O);!$(NS(=O)=O)]",
    "phys_amine_secondary": "[NX3;H1;!$(NC=O);!$(NS(=O)=O)]",
    "phys_amine_tertiary": "[NX3;H0;!$(NC=O);!$(NS(=O)=O)]",
    "phys_aromatic_n": "[nH0]",
    "phys_cationic_n": "[N+]",
}


def descriptor_row(smiles: str, patterns: dict[str, Chem.Mol]) -> dict[str, float]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid standardized SMILES: {smiles}")
    atoms = list(mol.GetAtoms())
    row = {
        "phys_mw": Descriptors.MolWt(mol),
        "phys_logp": Descriptors.MolLogP(mol),
        "phys_tpsa": Descriptors.TPSA(mol),
        "phys_hbd": rdMolDescriptors.CalcNumHBD(mol),
        "phys_hba": rdMolDescriptors.CalcNumHBA(mol),
        "phys_rotb": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "phys_rings": rdMolDescriptors.CalcNumRings(mol),
        "phys_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "phys_aliphatic_rings": rdMolDescriptors.CalcNumAliphaticRings(mol),
        "phys_frac_sp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "phys_heavy_atoms": mol.GetNumHeavyAtoms(),
        "phys_formal_charge": Chem.GetFormalCharge(mol),
        "phys_n_count": sum(atom.GetAtomicNum() == 7 for atom in atoms),
        "phys_o_count": sum(atom.GetAtomicNum() == 8 for atom in atoms),
        "phys_s_count": sum(atom.GetAtomicNum() == 16 for atom in atoms),
        "phys_halogen_count": sum(atom.GetAtomicNum() in (9, 17, 35, 53) for atom in atoms),
        "phys_amide": Fragments.fr_amide(mol),
        "phys_carboxylic_acid": Fragments.fr_COO(mol),
        "phys_phenol": Fragments.fr_phenol(mol),
        "phys_ether": Fragments.fr_ether(mol),
        "phys_ester": Fragments.fr_ester(mol),
        "phys_nitro": Fragments.fr_nitro(mol),
    }
    row.update({name: len(mol.GetSubstructMatches(pattern)) for name, pattern in patterns.items()})
    return row


def st_rae(frame: pl.DataFrame) -> float:
    y = frame["y_true"].to_numpy()
    pred = frame["prediction"].to_numpy()
    lo = frame["conf_low"].to_numpy()
    hi = frame["conf_high"].to_numpy()
    center = y.mean()
    return float(
        (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
        / (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    )


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    patterns = {name: Chem.MolFromSmarts(smarts) for name, smarts in AMINE_SMARTS.items()}
    if any(pattern is None for pattern in patterns.values()):
        raise ValueError("Invalid amine SMARTS")
    compounds = pl.read_parquet(inp / "compounds.parquet").sort("compound_id")
    physchem = pl.DataFrame(
        [
            {"compound_id": compound_id, **descriptor_row(smiles, patterns)}
            for compound_id, smiles in compounds.iter_rows()
        ]
    )
    phys_cols = sorted(c for c in physchem.columns if c.startswith("phys_"))
    if not np.isfinite(physchem.select(phys_cols).to_numpy()).all():
        raise ValueError("Nonfinite physchem features")
    physchem.write_parquet(out / "physchem.parquet")
    scores = pl.read_parquet(inp / "scores.parquet")
    medium = pl.read_parquet(inp / "medium.parquet")
    score_cols = sorted(c for c in scores.columns if c.startswith("score_"))
    medium_cols = sorted(c for c in medium.columns if c.startswith("emb_"))
    if len(score_cols) != 47 or len(medium_cols) != 600:
        raise ValueError("Unexpected frozen feature dimensions")
    features = physchem.join(
        scores.select("compound_id", *score_cols), on="compound_id", validate="1:1"
    ).join(medium.select("compound_id", *medium_cols), on="compound_id", validate="1:1")
    views = {
        "frozen": score_cols + medium_cols,
        "frozen_physchem": score_cols + medium_cols + phys_cols,
        "scores_physchem": score_cols + phys_cols,
        "physchem_only": phys_cols,
    }
    labels = (
        pl.read_parquet(inp / "measurements.parquet")
        .filter(pl.col("endpoint") == cfg["endpoint"])
        .select("compound_id", "endpoint", pl.col("value").alias("y_true"), "conf_low", "conf_high")
    )
    if labels.height != 1493 or list(views) != cfg["views"]:
        raise ValueError("Unexpected D6 labels or feature views")
    outputs = []
    for split in cfg["splits"]:
        assignments = pl.read_parquet(inp / f"{split}-assignments.parquet").select(
            "compound_id", "fold", "group", "is_anchor"
        )
        data = (
            labels.join(assignments, on="compound_id", validate="1:1")
            .join(features, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        if data.height != 1493 or data.filter(~pl.col("is_anchor")).height != 1457:
            raise ValueError(f"Incomplete {split} labels")
        fold_ids = data["fold"].to_numpy()
        anchors = data["is_anchor"].to_numpy()
        y = data["y_true"].to_numpy()
        for fold in range(5):
            train = np.flatnonzero((fold_ids != fold) | anchors)
            held = np.flatnonzero((fold_ids == fold) & ~anchors)
            for view in cfg["views"]:
                x = data.select(views[view]).to_numpy().astype(np.float32)
                if not np.isfinite(x).all():
                    raise ValueError(f"Nonfinite {split}/{view} features")
                params = cfg["model"]
                model = CatBoostRegressor(
                    iterations=params["iterations"],
                    depth=params["depth"],
                    learning_rate=params["learning_rate"],
                    l2_leaf_reg=params["l2_leaf_reg"],
                    loss_function=params["loss_function"],
                    thread_count=params["thread_count"],
                    random_seed=params["seed"] + fold,
                    verbose=False,
                    allow_writing_files=False,
                )
                model.fit(x[train], y[train])
                pred = model.predict(x[held])
                outputs.append(
                    data[held]
                    .select(
                        "compound_id",
                        "endpoint",
                        "fold",
                        "group",
                        "y_true",
                        "conf_low",
                        "conf_high",
                    )
                    .with_columns(
                        pl.lit(split).alias("split_id"),
                        pl.lit(view).alias("view"),
                        pl.Series("prediction", pred),
                    )
                )
                print(split, fold, view, flush=True)
    oof = pl.concat(outputs).sort("split_id", "view", "compound_id")
    if (
        oof.height != 2 * 4 * 1457
        or oof.select("split_id", "view", "compound_id").is_duplicated().any()
    ):
        raise ValueError("Incomplete OOF")
    metrics = []
    for (split, view), frame in oof.group_by("split_id", "view"):
        y = frame["y_true"].to_numpy()
        pred = frame["prediction"].to_numpy()
        metrics.append(
            {
                "split_id": split,
                "view": view,
                "n": frame.height,
                "ST_RAE": st_rae(frame),
                "MAE": float(np.abs(y - pred).mean()),
                "Spearman": float(spearmanr(y, pred).statistic),
                "mean_prediction": float(pred.mean()),
            }
        )
    oof.write_parquet(out / "oof.parquet")
    pl.DataFrame(metrics).sort("split_id", "view").write_csv(out / "metrics.csv")
    (out / "summary.json").write_text(
        json.dumps({"config": cfg, "physchem_columns": phys_cols}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
