"""Generate one-step CYP products and aggregate reaction/product features."""

import hashlib
import json
import os
import subprocess
import tarfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
scratch = Path(os.environ["TMPDIR"]) / "biotransformer-features"
scratch.mkdir(parents=True, exist_ok=False)
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
archive = inp / cfg["archive"]
if hashlib.sha256(archive.read_bytes()).hexdigest() != cfg["archive_sha256"]:
    raise ValueError("BioTransformer source archive hash mismatch")
with tarfile.open(archive) as handle:
    handle.extractall(scratch, filter="data")
jar = scratch / cfg["jar"]
if hashlib.sha256(jar.read_bytes()).hexdigest() != cfg["jar_sha256"]:
    raise ValueError("BioTransformer JAR hash mismatch")

compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
if cfg["max_compounds"] is not None:
    compounds = compounds.head(cfg["max_compounds"])
parent_by_key = {}
records = list(compounds.iter_rows(named=True))
for record in records:
    molecule = Chem.MolFromSmiles(record["standardized_smiles"])
    parent_by_key[Chem.MolToInchiKey(molecule)] = (record["compound_id"], molecule)


def run_chunk(chunk_index, chunk_records):
    input_sdf = scratch / f"challenge-{chunk_index:04d}.sdf"
    raw_output = scratch / f"products-{chunk_index:04d}.sdf"
    writer = Chem.SDWriter(str(input_sdf))
    for record in chunk_records:
        molecule = Chem.MolFromSmiles(record["standardized_smiles"])
        molecule.SetProp("_Name", record["compound_id"])
        molecule.SetProp("compound_id", record["compound_id"])
        writer.write(molecule)
    writer.close()
    command = [
        "java",
        f"-Xmx{cfg['java_heap']}",
        "-jar",
        str(jar),
        "-k",
        "pred",
        "-b",
        "cyp450",
        "-isdf",
        str(input_sdf),
        "-osdf",
        str(raw_output),
        "-s",
        str(cfg["steps"]),
        "-cm",
        str(cfg["cyp450_mode"]),
    ]
    completed = subprocess.run(command, cwd=scratch, text=True, capture_output=True, check=True)
    return chunk_index, raw_output, command, completed.stdout, completed.stderr


chunks = [
    records[start : start + cfg["chunk_size"]]
    for start in range(0, len(records), cfg["chunk_size"])
]
results = []
with ThreadPoolExecutor(max_workers=cfg["workers"]) as executor:
    futures = [executor.submit(run_chunk, index, chunk) for index, chunk in enumerate(chunks)]
    for future in as_completed(futures):
        result = future.result()
        results.append(result)
        print(f"BioTransformer chunk {result[0] + 1}/{len(chunks)} complete", flush=True)
results.sort()
(out / "stdout.txt").write_text("\n".join(result[3] for result in results))
(out / "stderr.txt").write_text("\n".join(result[4] for result in results))
raw_output = out / "products.sdf"
raw_writer = Chem.SDWriter(str(raw_output))
for _, chunk_output, _, _, _ in results:
    for molecule in Chem.SDMolSupplier(str(chunk_output), removeHs=False):
        if molecule is not None:
            raw_writer.write(molecule)
raw_writer.close()
commands = [result[2] for result in results]

products_by_id = defaultdict(list)
unmapped = []
for product in Chem.SDMolSupplier(str(raw_output), removeHs=False):
    if product is None:
        unmapped.append("invalid_product")
        continue
    key = product.GetProp("Precursor_InChIKey") if product.HasProp("Precursor_InChIKey") else ""
    if key not in parent_by_key:
        precursor = (
            Chem.MolFromSmiles(product.GetProp("Precursor_SMILES"))
            if product.HasProp("Precursor_SMILES")
            else None
        )
        key = Chem.MolToInchiKey(precursor) if precursor is not None else ""
    if key not in parent_by_key:
        unmapped.append(key)
        continue
    products_by_id[parent_by_key[key][0]].append(product)

keywords = ["hydroxyl", "dealkyl", "desatur", "epoxid", "oxid", "reduct", "cleavage"]
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=cfg["fingerprint_bits"])
feature_rows, fingerprint_rows, product_rows = [], [], []
for record in compounds.iter_rows(named=True):
    compound_id = record["compound_id"]
    parent = Chem.MolFromSmiles(record["standardized_smiles"])
    products = products_by_id[compound_id]
    unique = {}
    for product in products:
        unique.setdefault(Chem.MolToSmiles(product), product)
    product_list = list(unique.values())
    enzymes = [
        p.GetProp("Enzyme(s)").splitlines() if p.HasProp("Enzyme(s)") else [] for p in product_list
    ]
    reactions = [
        p.GetProp("Reaction").lower() if p.HasProp("Reaction") else "" for p in product_list
    ]
    masses = np.asarray([Descriptors.ExactMolWt(p) for p in product_list], dtype=float)
    parent_mass = Descriptors.ExactMolWt(parent)
    parent_fp = generator.GetFingerprintAsNumPy(parent).astype(np.float32)
    product_fps = (
        np.stack([generator.GetFingerprintAsNumPy(p) for p in product_list]).astype(np.float32)
        if product_list
        else np.zeros((0, cfg["fingerprint_bits"]), dtype=np.float32)
    )
    mean_fp = product_fps.mean(axis=0) if len(product_fps) else np.zeros_like(parent_fp)
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
            **{f"biotransformer_product_fp_{i}": float(value) for i, value in enumerate(mean_fp)},
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
products_frame = (
    pl.DataFrame(product_rows)
    if product_rows
    else pl.DataFrame(
        schema={
            "compound_id": pl.String,
            "product_smiles": pl.String,
            "reaction": pl.String,
            "reaction_id": pl.String,
            "enzymes": pl.String,
        }
    )
)
products_frame.write_parquet(out / "biotransformer-products.parquet")
summary = {
    "source_commit": cfg["source_commit"],
    "commands": commands,
    "attempted_compounds": compounds.height,
    "compounds_with_products": len(products_by_id),
    "unique_products": products_frame.height,
    "unmapped_products": len(unmapped),
    "product_fingerprint_pooling": (
        "mean Morgan radius-2 bit vector over unique predicted metabolites"
    ),
    "workers": cfg["workers"],
    "chunk_size": cfg["chunk_size"],
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
