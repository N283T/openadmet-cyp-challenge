"""Extract frozen parent and mean BioTransformer-product CheMeleon embeddings."""

import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop import nn as cnn
from chemprop.data import BatchMolGraph
from rdkit import Chem

from cyp_challenge.storage import digest

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
scratch = Path(os.environ["TMPDIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
products = pl.read_parquet(inp / cfg["products"]).select(
    "compound_id", "product_smiles", "enzymes"
).unique()
checkpoint_path = inp / cfg["checkpoint"]
if digest(checkpoint_path) != cfg["checkpoint_sha256"]:
    raise ValueError("CheMeleon checkpoint hash mismatch")
checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
torch.manual_seed(cfg["seed"])
mp = cnn.BondMessagePassing(**checkpoint["hyper_parameters"])
mp.load_state_dict(checkpoint["state_dict"], strict=True)
mp.to(device).eval()
aggregation = cnn.MeanAggregation()
featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()

parent_smiles = dict(compounds.select("compound_id", "standardized_smiles").iter_rows())
structure_smiles = sorted(
    set(compounds["standardized_smiles"].to_list()) | set(products["product_smiles"].to_list())
)
smiles_to_index = {smiles: index for index, smiles in enumerate(structure_smiles)}
dimension = mp.output_dim
embedding_path = scratch / "structure-embeddings.float32"
embeddings = np.memmap(
    embedding_path, mode="w+", dtype=np.float32, shape=(len(structure_smiles), dimension)
)

with torch.inference_mode():
    for start in range(0, len(structure_smiles), cfg["batch_size"]):
        batch_smiles = structure_smiles[start : start + cfg["batch_size"]]
        graphs = [featurizer(Chem.MolFromSmiles(smiles)) for smiles in batch_smiles]
        graph = BatchMolGraph(graphs)
        graph.to(device)
        values = aggregation(mp(graph), graph.batch).cpu().numpy().astype(np.float32)
        embeddings[start : start + len(batch_smiles)] = values
        if start % (cfg["batch_size"] * 50) == 0:
            print(f"embedded {start + len(batch_smiles)}/{len(structure_smiles)}", flush=True)
embeddings.flush()

scope_indices: dict[str, dict[str, list[int]]] = defaultdict(
    lambda: {"all": [], "d6": [], "a4": []}
)
for record in products.iter_rows(named=True):
    index = smiles_to_index[record["product_smiles"]]
    scope_indices[record["compound_id"]]["all"].append(index)
    enzymes = record["enzymes"] or ""
    if "CYP2D6" in enzymes:
        scope_indices[record["compound_id"]]["d6"].append(index)
    if "CYP3A4" in enzymes:
        scope_indices[record["compound_id"]]["a4"].append(index)

n_compounds = compounds.height
blocks = {
    scope: np.empty((n_compounds, dimension), dtype=np.float32)
    for scope in ["parent", "all", "d6", "a4"]
}
scalar_rows = []
for row, compound_id in enumerate(compounds["compound_id"]):
    parent = np.asarray(embeddings[smiles_to_index[parent_smiles[compound_id]]])
    blocks["parent"][row] = parent
    scalar = {"compound_id": compound_id}
    parent_norm = max(float(np.linalg.norm(parent)), 1e-12)
    for scope in ["all", "d6", "a4"]:
        indices = sorted(set(scope_indices[compound_id][scope]))
        mean = np.asarray(embeddings[indices]).mean(axis=0) if indices else parent
        blocks[scope][row] = mean
        mean_norm = max(float(np.linalg.norm(mean)), 1e-12)
        scalar[f"chemeleon_{scope}_count"] = len(indices)
        scalar[f"chemeleon_{scope}_parent_cosine"] = float(
            np.dot(parent, mean) / (parent_norm * mean_norm)
        )
        scalar[f"chemeleon_{scope}_parent_l2"] = float(np.linalg.norm(mean - parent))
    scalar_rows.append(scalar)

columns = {"compound_id": compounds["compound_id"]}
for scope, matrix in blocks.items():
    for index in range(dimension):
        columns[f"chemeleon_{scope}_{index:04d}"] = matrix[:, index]
pl.DataFrame(columns).write_parquet(out / "product-embeddings.parquet")
pl.DataFrame(scalar_rows).write_parquet(out / "product-embedding-scalars.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "selection": cfg["selection"],
            "source_run": cfg["source_run"],
            "device": str(device),
            "structures_embedded": len(structure_smiles),
            "parent_product_pairs": products.height,
            "dimension": dimension,
            "blocks": list(blocks),
        },
        indent=2,
    )
    + "\n"
)
print((out / "summary.json").read_text(), flush=True)
