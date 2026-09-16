"""Score BioTransformer products with fold-safe direct and log2fc regressors."""

import json
import os
from collections import defaultdict
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


compounds = read("compounds").sort("compound_id")
measurements = read("measurements")
targets = read("log2fc_targets")
assignments = read("assignments")
products = read("products").select("compound_id", "product_smiles", "enzymes").unique()
test_ids = set(read("test_mapping")["compound_id"])

parent_smiles = dict(compounds.select("compound_id", "standardized_smiles").iter_rows())
structure_smiles = sorted(
    set(compounds["standardized_smiles"].to_list()) | set(products["product_smiles"].to_list())
)
smiles_to_index = {smiles: index for index, smiles in enumerate(structure_smiles)}
generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=cfg["fingerprint_radius"], fpSize=cfg["fingerprint_bits"]
)
fingerprints = np.stack(
    [generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles)) for smiles in structure_smiles]
).astype(np.uint8)
parent_indices = {
    compound_id: smiles_to_index[smiles] for compound_id, smiles in parent_smiles.items()
}

product_indices: dict[str, dict[str, list[int]]] = defaultdict(
    lambda: {"all": [], "d6": [], "a4": []}
)
for record in products.iter_rows(named=True):
    index = smiles_to_index[record["product_smiles"]]
    product_indices[record["compound_id"]]["all"].append(index)
    enzymes = record["enzymes"] or ""
    if "CYP2D6" in enzymes:
        product_indices[record["compound_id"]]["d6"].append(index)
    if "CYP3A4" in enzymes:
        product_indices[record["compound_id"]]["a4"].append(index)
for compound_id in parent_smiles:
    for scope in ["all", "d6", "a4"]:
        product_indices[compound_id][scope] = sorted(set(product_indices[compound_id][scope]))

auxiliary_targets = {}
for enzyme in cfg["enzymes"]:
    direct_endpoint = f"{enzyme}_pIC50_direct_inhibition"
    direct = measurements.filter(pl.col("endpoint") == direct_endpoint).select(
        "compound_id", "value"
    )
    log_head = next(head for head in targets["head"].unique() if head.startswith(f"{enzyme}@"))
    log2fc = targets.filter(pl.col("head") == log_head).select("compound_id", "value")
    auxiliary_targets[f"{enzyme.lower()}_direct"] = direct
    auxiliary_targets[f"{enzyme.lower()}_log2fc"] = log2fc


def aggregate(scores: np.ndarray, compound_id: str, prefix: str) -> dict[str, float]:
    parent = float(scores[parent_indices[compound_id]])
    row = {f"{prefix}_parent": parent}
    for scope in ["all", "d6", "a4"]:
        indices = product_indices[compound_id][scope]
        values = scores[indices] if indices else np.asarray([parent])
        ordered = np.sort(values)[::-1]
        row[f"{prefix}_{scope}_mean"] = float(values.mean())
        row[f"{prefix}_{scope}_std"] = float(values.std())
        row[f"{prefix}_{scope}_min"] = float(values.min())
        row[f"{prefix}_{scope}_max"] = float(values.max())
        row[f"{prefix}_{scope}_top3_mean"] = float(ordered[:3].mean())
        row[f"{prefix}_{scope}_max_delta"] = float(values.max() - parent)
        row[f"{prefix}_{scope}_min_delta"] = float(values.min() - parent)
        row[f"{prefix}_{scope}_mean_delta"] = float(values.mean() - parent)
    return row


feature_rows = []
fit_rows = []
for strategy in cfg["strategies"]:
    fold_table = assignments.filter(pl.col("strategy") == strategy).select(
        "compound_id", "fold"
    )
    for fold in range(5):
        score_by_target = {}
        for target_name, label_table in auxiliary_targets.items():
            training = (
                label_table.join(fold_table, on="compound_id", validate="1:1")
                .filter(pl.col("fold") != fold)
                .sort("compound_id")
            )
            train_indices = np.asarray([parent_indices[cid] for cid in training["compound_id"]])
            model = lgb.LGBMRegressor(**cfg["model"], random_state=cfg["seed"] + fold)
            model.fit(fingerprints[train_indices], training["value"].to_numpy())
            score_by_target[target_name] = model.predict(fingerprints).astype(np.float32)
            model.booster_.save_model(str(out / f"{strategy}-{fold}-{target_name}.lgbm.txt"))
            fit_rows.append(
                {
                    "strategy": strategy,
                    "fold": fold,
                    "target": target_name,
                    "training_rows": training.height,
                }
            )
        for compound_id in compounds["compound_id"]:
            row = {"compound_id": compound_id, "strategy": strategy, "model_fold": fold}
            for target_name, scores in score_by_target.items():
                row.update(aggregate(scores, compound_id, target_name))
            feature_rows.append(row)
        print(f"{strategy} fold {fold} complete", flush=True)

features = pl.DataFrame(feature_rows).sort("strategy", "model_fold", "compound_id")
if features.height != len(cfg["strategies"]) * 5 * compounds.height:
    raise ValueError("Incomplete strategy/fold/compound feature grid")
if set(features.filter(pl.col("compound_id").is_in(test_ids))["compound_id"]) != test_ids:
    raise ValueError("Test feature coverage is incomplete")
features.write_parquet(out / "product-activity-features.parquet")
pl.DataFrame(fit_rows).write_csv(out / "fits.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "selection": cfg["selection"],
            "structures_scored": len(structure_smiles),
            "parent_product_pairs": products.height,
            "feature_rows": features.height,
            "feature_count": len(features.columns) - 3,
            "fit_count": len(fit_rows),
            "routing": "Downstream outer fold uses the auxiliary transform excluding that fold",
        },
        indent=2,
    )
    + "\n"
)
print((out / "summary.json").read_text(), flush=True)
