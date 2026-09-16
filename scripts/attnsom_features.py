"""Extract ATTNSOM SoM probabilities and pretrained molecular representations."""

import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np
import polars as pl
import torch
from rdkit import Chem
from torch_geometric.loader import DataLoader

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


dataset_utils = load_module("dataset_utils", inp / cfg["dataset_utils"])
model_module = load_module("attnsom_model", inp / cfg["model_source"])
model_config = json.loads((inp / cfg["model_config"]).read_text())
device = torch.device(cfg["device"] if torch.cuda.is_available() else "cpu")
cyp_index = {name: index for index, name in enumerate(model_config["cyp_names"])}

compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
graphs = []
graph_ids = []
failures = []
for record in compounds.iter_rows(named=True):
    molecule = Chem.MolFromSmiles(record["standardized_smiles"])
    if molecule is None or molecule.GetNumAtoms() == 0:
        failures.append({"compound_id": record["compound_id"], "reason": "rdkit_parse"})
        continue
    try:
        graph, _ = dataset_utils.mol_to_graph(molecule, Chem.MolToSmiles(molecule), som_indices=[])
    except Exception as error:
        failures.append(
            {"compound_id": record["compound_id"], "reason": f"{type(error).__name__}: {error}"}
        )
        continue
    graph.compound_position = torch.tensor([len(graphs)], dtype=torch.long)
    graphs.append(graph)
    graph_ids.append(record["compound_id"])

if not graphs:
    raise ValueError("No challenge molecules could be featurized")

n_molecules = len(graphs)
n_folds = len(cfg["checkpoints"])
embedding_sum = np.zeros((n_molecules, 2 * model_config["hidden_size"]), dtype=np.float64)
probability_sum = {
    isoform: [np.zeros(graph.num_nodes, dtype=np.float64) for graph in graphs]
    for isoform in cfg["isoforms"]
}

for fold_index, checkpoint_name in enumerate(cfg["checkpoints"], start=1):
    model = model_module.ATTNSOM(**model_config).to(device)
    checkpoint = torch.load(inp / checkpoint_name, map_location=device, weights_only=False)
    state = checkpoint.get("model_state_dict", checkpoint.get("state_dict", checkpoint))
    model.load_state_dict(state)
    model.eval()
    for isoform_index, isoform in enumerate(cfg["isoforms"]):
        for graph in graphs:
            graph.cyp_idx = torch.full((graph.num_nodes,), cyp_index[isoform], dtype=torch.long)
        loader = DataLoader(graphs, batch_size=cfg["batch_size"], shuffle=False)
        with torch.inference_mode():
            for batch in loader:
                batch = batch.to(device)
                logits, h_modulated, _ = model(batch)
                probabilities = torch.sigmoid(logits).cpu().numpy()
                representations = h_modulated.cpu().numpy()
                batch_index = batch.batch.cpu().numpy()
                positions = batch.compound_position.cpu().numpy()
                offset = 0
                for local_index, position in enumerate(positions):
                    mask = batch_index == local_index
                    count = int(mask.sum())
                    probability_sum[isoform][position] += probabilities[offset : offset + count]
                    if isoform_index == 0:
                        atom_repr = representations[offset : offset + count]
                        embedding_sum[position] += np.concatenate(
                            [atom_repr.mean(axis=0), atom_repr.max(axis=0)]
                        )
                    offset += count
    print(f"ATTNSOM extracted fold {fold_index}/{n_folds}", flush=True)

embeddings = embedding_sum / n_folds
embedding_rows = []
feature_rows = []
atom_rows = []
for position, (compound_id, graph) in enumerate(zip(graph_ids, graphs, strict=True)):
    embedding_rows.append(
        {
            "compound_id": compound_id,
            **{
                f"attnsom_emb_{index}": float(value)
                for index, value in enumerate(embeddings[position])
            },
        }
    )
    feature_row = {"compound_id": compound_id}
    for isoform in cfg["isoforms"]:
        values = probability_sum[isoform][position] / n_folds
        ordered = np.sort(values)[::-1]
        prefix = f"attnsom_{isoform.lower()}"
        feature_row.update(
            {
                f"{prefix}_max": float(ordered[0]),
                f"{prefix}_top2_mean": float(ordered[:2].mean()),
                f"{prefix}_top3_mean": float(ordered[:3].mean()),
                f"{prefix}_mean": float(values.mean()),
                f"{prefix}_max_gap": float(ordered[0] - ordered[1]) if len(ordered) > 1 else 0.0,
                f"{prefix}_count_gt_025": int((values > 0.25).sum()),
                f"{prefix}_count_gt_050": int((values > 0.50).sum()),
            }
        )
        molecule = Chem.MolFromSmiles(graph.smiles)
        for atom_index, probability in enumerate(values):
            atom = molecule.GetAtomWithIdx(atom_index)
            atom_rows.append(
                {
                    "compound_id": compound_id,
                    "isoform": isoform,
                    "atom_index": atom_index,
                    "atomic_number": atom.GetAtomicNum(),
                    "probability": float(probability),
                    "rank": int((-values).argsort().argsort()[atom_index] + 1),
                }
            )
    feature_rows.append(feature_row)

pl.DataFrame(embedding_rows).write_parquet(out / "attnsom-embeddings.parquet")
pl.DataFrame(feature_rows).write_parquet(out / "attnsom-features.parquet")
pl.DataFrame(atom_rows).write_parquet(out / "attnsom-atoms.parquet")
failure_frame = (
    pl.DataFrame(failures)
    if failures
    else pl.DataFrame(schema={"compound_id": pl.String, "reason": pl.String})
)
failure_frame.write_parquet(out / "attnsom-failures.parquet")
summary = {
    "source_commit": cfg["source_commit"],
    "training_run": cfg["training_run"],
    "training_checkpoints": n_folds,
    "attempted_compounds": compounds.height,
    "successful_compounds": n_molecules,
    "failed_compounds": failure_frame.height,
    "isoforms": cfg["isoforms"],
    "embedding_width": embeddings.shape[1],
    "pooling": (
        "fold-mean of atomwise mean and componentwise maximum of "
        "FiLM-modulated representations"
    ),
    "leakage_control": "all challenge molecular connectivities were excluded from ATTNSOM training",
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
