"""Extract frozen CYP-MAP substrate, SoM, reaction-type, and embedding features."""

import json
import os
import re
import sys
import types
from pathlib import Path

import numpy as np
import polars as pl
import torch
from torch_geometric.data import Batch, Data
from torch_geometric.utils import scatter as pyg_scatter

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
source = Path(cfg["source_directory"]).expanduser()
sys.path.insert(0, str(source))

# CYP-MAP imports torch_scatter only for this operation; modern PyG provides the
# same scatter API and avoids pinning an obsolete compiled extension.
scatter_module = types.ModuleType("torch_scatter")
scatter_module.scatter = pyg_scatter
sys.modules.setdefault("torch_scatter", scatter_module)

from modules.dualgraph.mol import smiles2graphwithface  # noqa: E402
from modules.som_models import CYPMAP_GNN  # noqa: E402


class CYPMapData(Data):
    """CYP-MAP graph with the upstream batching offsets."""

    def __cat_dim__(self, key, value, *args, **kwargs):
        if bool(re.search("(index|face)", key)) or bool(re.search("(nf_node|nf_ring)", key)):
            return -1
        return 0

    def __inc__(self, key, value, *args, **kwargs):
        if bool(re.search("(ring_index|nf_ring)", key)):
            return int(self.num_rings)
        if bool(re.search("(index|face|nf_node)", key)):
            return int(self.num_nodes)
        return 0


def make_graph(smiles: str) -> CYPMapData:
    graph = smiles2graphwithface(smiles)
    data = CYPMapData()
    data.num_nodes = int(graph["num_nodes"])
    data.x = torch.from_numpy(graph["node_feat"]).long()
    data.edge_index = torch.from_numpy(graph["edge_index"]).long()
    data.edge_attr = torch.from_numpy(graph["edge_feat"]).long()
    data.ring_mask = torch.from_numpy(graph["ring_mask"]).bool()
    data.ring_index = torch.from_numpy(graph["ring_index"]).long()
    data.nf_node = torch.from_numpy(graph["nf_node"]).long()
    data.nf_ring = torch.from_numpy(graph["nf_ring"]).long()
    data.num_rings = int(graph["num_rings"])
    data.n_edges = int(graph["n_edges"])
    data.n_nodes = int(graph["n_nodes"])
    data.n_nfs = int(graph["n_nfs"])
    return data


def summarize(values: torch.Tensor) -> tuple[float, float, float]:
    if values.numel() == 0:
        return 0.0, 0.0, 0.0
    values = values.sigmoid()
    k = min(3, values.numel())
    return float(values.max()), float(values.mean()), float(values.topk(k).values.mean())


compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
if cfg.get("max_molecules"):
    compounds = compounds.head(cfg["max_molecules"])

graphs = []
valid_ids = []
failed = []
for compound_id, smiles in compounds.select("compound_id", "standardized_smiles").iter_rows():
    try:
        graph = make_graph(smiles)
        if graph.num_nodes == 0 or graph.edge_index.shape[1] == 0:
            raise ValueError("CYP-MAP requires at least one atom and one bond")
        graphs.append(graph)
        valid_ids.append(compound_id)
    except Exception as exc:
        failed.append({"compound_id": compound_id, "error": f"{type(exc).__name__}: {exc}"})

device = torch.device("cuda" if cfg["device"] == "cuda" and torch.cuda.is_available() else "cpu")
cyp_names = {"2D6": "BOM_2D6", "3A4": "BOM_3A4"}
site_tasks = {
    "atom": [
        "atom_som",
        "atom_spn",
        "atom_oxc",
        "atom_oxi",
        "atom_epo",
        "atom_sut",
        "atom_dhy",
        "atom_hys",
        "atom_rdc",
    ],
    "bond": [
        "bond_som",
        "bond_oxc",
        "bond_oxi",
        "bond_epo",
        "bond_sut",
        "bond_dhy",
        "bond_hys",
        "bond_rdc",
    ],
}
model_feature_rows = []
model_embeddings = []

for model_index, checkpoint_name in enumerate(cfg["checkpoints"]):
    model = CYPMAP_GNN(
        cyp_list=[
            "BOM_1A2",
            "BOM_2A6",
            "BOM_2B6",
            "BOM_2C8",
            "BOM_2C9",
            "BOM_2C19",
            "BOM_2D6",
            "BOM_2E1",
            "BOM_3A4",
            "CYP_REACTION",
        ]
    ).to(device)
    state = torch.load(inp / checkpoint_name, map_location="cpu", weights_only=True)
    model.load_state_dict(state)
    model.eval()
    captured = []

    def capture_embedding(module, args, sink=captured):
        sink.append(args[0].detach().cpu())

    handle = model.substrate_fc.register_forward_pre_hook(capture_embedding)
    rows = []
    embeddings = []
    with torch.inference_mode():
        for start in range(0, len(graphs), cfg["batch_size"]):
            items = graphs[start : start + cfg["batch_size"]]
            batch = Batch.from_data_list(items).to(device)
            captured.clear()
            logits = model(batch)
            molecular = captured[0]
            atom_counts = [item.num_nodes for item in items]
            bond_counts = [item.edge_index.shape[1] // 2 for item in items]
            atom_offsets = np.cumsum([0, *atom_counts])
            bond_offsets = np.cumsum([0, *bond_counts])
            for local_index, compound_id in enumerate(valid_ids[start : start + len(items)]):
                row = {"compound_id": compound_id, "model_index": model_index}
                for short_cyp, cyp in cyp_names.items():
                    prefix = short_cyp.lower()
                    row[f"{prefix}_substrate"] = float(
                        logits["subs"][cyp][local_index].sigmoid().cpu()
                    )
                    for level, tasks in site_tasks.items():
                        lo, hi = (
                            (atom_offsets[local_index], atom_offsets[local_index + 1])
                            if level == "atom"
                            else (bond_offsets[local_index], bond_offsets[local_index + 1])
                        )
                        for task in tasks:
                            maximum, mean, top3 = summarize(logits[task][cyp][lo:hi].cpu())
                            name = task.removeprefix(f"{level}_")
                            row[f"{prefix}_{level}_{name}_max"] = maximum
                            if name == "som":
                                row[f"{prefix}_{level}_{name}_mean"] = mean
                                row[f"{prefix}_{level}_{name}_top3"] = top3
                rows.append(row)
                embeddings.append(molecular[local_index].numpy().astype(np.float32))
            progress = min(start + len(items), len(graphs))
            print(
                f"model {model_index + 1}/{len(cfg['checkpoints'])}: "
                f"{progress}/{len(graphs)}",
                flush=True,
            )
    handle.remove()
    model_feature_rows.append(pl.DataFrame(rows).sort("compound_id"))
    model_embeddings.append(np.stack(embeddings))
    del model, state
    if device.type == "cuda":
        torch.cuda.empty_cache()

id_table = pl.DataFrame({"compound_id": valid_ids}).sort("compound_id")
numeric_columns = [
    column
    for column in model_feature_rows[0].columns
    if column not in {"compound_id", "model_index"}
]
feature_stack = np.stack(
    [table.select(numeric_columns).to_numpy().astype(np.float32) for table in model_feature_rows]
)
feature_data = {"compound_id": id_table["compound_id"]}
for column_index, column in enumerate(numeric_columns):
    feature_data[f"cypmap_{column}_mean"] = feature_stack[:, :, column_index].mean(axis=0)
    feature_data[f"cypmap_{column}_std"] = feature_stack[:, :, column_index].std(axis=0)
features = pl.DataFrame(feature_data)
features.write_parquet(out / "cypmap-features.parquet")

embedding_stack = np.stack(model_embeddings)
embedding_mean = embedding_stack.mean(axis=0)
embedding_data = {"compound_id": id_table["compound_id"]}
for index in range(embedding_mean.shape[1]):
    embedding_data[f"cypmap_emb_{index:03d}"] = embedding_mean[:, index]
pl.DataFrame(embedding_data).write_parquet(out / "cypmap-embeddings.parquet")

failure_table = (
    pl.DataFrame(failed)
    if failed
    else pl.DataFrame(schema={"compound_id": pl.String, "error": pl.String})
)
failure_table.write_parquet(out / "cypmap-failures.parquet")
summary = {
    "source_commit": cfg["source_commit"],
    "device": str(device),
    "models": len(cfg["checkpoints"]),
    "requested_molecules": compounds.height,
    "successful_molecules": len(valid_ids),
    "failed_molecules": len(failed),
    "property_features": len(features.columns) - 1,
    "embedding_dimensions": embedding_mean.shape[1],
    "representation": (
        "Mean across five frozen CYP-MAP trials; molecular embedding is the 512-dimensional "
        "input to the nine-CYP substrate head. Property features include per-CYP substrate, "
        "atom/bond SoM, reaction-type maxima, and ensemble standard deviations."
    ),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
