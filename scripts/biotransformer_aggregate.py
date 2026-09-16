"""Map BioTransformer products robustly and create molecular feature tables."""

import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")

by_smiles = {}
by_key = {}
by_connectivity = defaultdict(list)
parent_molecules = {}
for record in compounds.iter_rows(named=True):
    molecule = Chem.MolFromSmiles(record["standardized_smiles"])
    compound_id = record["compound_id"]
    canonical = Chem.MolToSmiles(molecule)
    key = Chem.MolToInchiKey(molecule)
    by_smiles[canonical] = compound_id
    by_key[key] = compound_id
    by_connectivity[key.split("-", 1)[0]].append(compound_id)
    parent_molecules[compound_id] = molecule

products_by_id = defaultdict(list)
unmapped_rows = []
mapping_counts = defaultdict(int)
for source_row, product in enumerate(
    Chem.SDMolSupplier(str(inp / cfg["products"]), removeHs=False)
):
    if product is None:
        unmapped_rows.append({"source_row": source_row, "reason": "invalid_product"})
        continue
    precursor_smiles = (
        product.GetProp("Precursor_SMILES") if product.HasProp("Precursor_SMILES") else ""
    )
    precursor = Chem.MolFromSmiles(precursor_smiles)
    canonical = Chem.MolToSmiles(precursor) if precursor is not None else ""
    supplied_key = (
        product.GetProp("Precursor_InChIKey") if product.HasProp("Precursor_InChIKey") else ""
    )
    compound_id = None
    method = None
    if canonical in by_smiles:
        compound_id, method = by_smiles[canonical], "canonical_smiles"
    elif supplied_key in by_key:
        compound_id, method = by_key[supplied_key], "full_inchikey"
    else:
        key = Chem.MolToInchiKey(precursor) if precursor is not None else supplied_key
        candidates = by_connectivity.get(key.split("-", 1)[0], []) if key else []
        if len(candidates) == 1:
            compound_id, method = candidates[0], "unique_connectivity"
    if compound_id is None:
        unmapped_rows.append(
            {
                "source_row": source_row,
                "reason": "ambiguous_or_absent_parent",
                "precursor_smiles": precursor_smiles,
                "precursor_inchikey": supplied_key,
                "product_smiles": Chem.MolToSmiles(product),
            }
        )
        continue
    mapping_counts[method] += 1
    products_by_id[compound_id].append(product)

keywords = ["hydroxyl", "dealkyl", "desatur", "epoxid", "oxid", "reduct", "cleavage"]
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=cfg["fingerprint_bits"])
feature_rows = []
fingerprint_rows = []
product_rows = []
for compound_id in compounds["compound_id"]:
    parent = parent_molecules[compound_id]
    unique = {}
    for product in products_by_id[compound_id]:
        unique.setdefault(Chem.MolToSmiles(product), product)
    product_list = list(unique.values())
    enzymes = [
        product.GetProp("Enzyme(s)").splitlines() if product.HasProp("Enzyme(s)") else []
        for product in product_list
    ]
    reactions = [
        product.GetProp("Reaction").lower() if product.HasProp("Reaction") else ""
        for product in product_list
    ]
    masses = np.asarray([Descriptors.ExactMolWt(product) for product in product_list])
    parent_mass = Descriptors.ExactMolWt(parent)
    product_fps = (
        np.stack([generator.GetFingerprintAsNumPy(product) for product in product_list])
        if product_list
        else np.zeros((0, cfg["fingerprint_bits"]), dtype=np.float32)
    )
    mean_fp = (
        product_fps.mean(axis=0)
        if len(product_fps)
        else np.zeros(cfg["fingerprint_bits"], dtype=np.float32)
    )
    feature_rows.append(
        {
            "compound_id": compound_id,
            "biotransformer_product_count": len(product_list),
            "biotransformer_d6_product_count": sum("CYP2D6" in value for value in enzymes),
            "biotransformer_a4_product_count": sum("CYP3A4" in value for value in enzymes),
            "biotransformer_d6_a4_product_count": sum(
                "CYP2D6" in value and "CYP3A4" in value for value in enzymes
            ),
            "biotransformer_reaction_count": len(set(reactions)),
            "biotransformer_mass_delta_mean": float((masses - parent_mass).mean())
            if len(masses)
            else 0.0,
            "biotransformer_mass_delta_min": float((masses - parent_mass).min())
            if len(masses)
            else 0.0,
            "biotransformer_mass_delta_max": float((masses - parent_mass).max())
            if len(masses)
            else 0.0,
            **{
                f"biotransformer_{keyword}_count": sum(
                    keyword in reaction for reaction in reactions
                )
                for keyword in keywords
            },
        }
    )
    fingerprint_rows.append(
        {
            "compound_id": compound_id,
            **{
                f"biotransformer_product_fp_{index}": float(value)
                for index, value in enumerate(mean_fp)
            },
        }
    )
    for product in product_list:
        product_rows.append(
            {
                "compound_id": compound_id,
                "product_smiles": Chem.MolToSmiles(product),
                "reaction": product.GetProp("Reaction") if product.HasProp("Reaction") else None,
                "reaction_id": product.GetProp("Reaction_ID")
                if product.HasProp("Reaction_ID")
                else None,
                "enzymes": product.GetProp("Enzyme(s)") if product.HasProp("Enzyme(s)") else None,
            }
        )

pl.DataFrame(feature_rows).write_parquet(out / "biotransformer-features.parquet")
pl.DataFrame(fingerprint_rows).write_parquet(out / "biotransformer-product-fingerprints.parquet")
pl.DataFrame(product_rows).write_parquet(out / "biotransformer-products.parquet")
unmapped = (
    pl.DataFrame(unmapped_rows)
    if unmapped_rows
    else pl.DataFrame(schema={"source_row": pl.Int64, "reason": pl.String})
)
unmapped.write_parquet(out / "unmapped-products.parquet")
summary = {
    "generation_run": cfg["generation_run"],
    "attempted_compounds": compounds.height,
    "compounds_with_products": len(products_by_id),
    "unique_products": len(product_rows),
    "unmapped_products": unmapped.height,
    "mapping_counts": dict(mapping_counts),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
