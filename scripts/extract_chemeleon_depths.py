"""Export pooled CheMeleon states after selected message-passing iterations."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers, nn
from chemprop.data import BatchMolGraph
from rdkit import Chem


def states(mp, graph, depths):
    """Use Chemprop's forward steps; apply its final atom projection at each depth."""
    graph = mp.graph_transform(graph)
    initial = mp.initialize(graph)
    hidden = mp.tau(initial)
    result = {}
    for depth in range(1, mp.depth + 1):
        if depth in depths:
            index = graph.edge_index[1].unsqueeze(1).expand(-1, hidden.shape[1])
            messages = torch.zeros(
                len(graph.V), hidden.shape[1], dtype=hidden.dtype, device=hidden.device
            ).scatter_reduce_(0, index, hidden, reduce="sum", include_self=False)
            result[depth] = mp.finalize(messages, graph.V, None)
        if depth < mp.depth:
            current = (hidden + hidden[graph.rev_edge_index]) / 2 if mp.undirected else hidden
            hidden = mp.update(mp.message(current, graph), initial)
    return result


def main():
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    torch.set_num_threads(cfg["threads"])
    checkpoint = torch.load(inp / cfg["checkpoint"], map_location="cpu", weights_only=True)
    mp = nn.BondMessagePassing(**checkpoint["hyper_parameters"])
    mp.load_state_dict(checkpoint["state_dict"], strict=True)
    mp = mp.to(cfg["device"]).eval().requires_grad_(False)
    aggregation = nn.MeanAggregation().to(cfg["device"]).eval()
    depths = cfg["depths"]
    assert depths == sorted(set(depths)) and depths[-1] == mp.depth
    compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
    assert compounds["compound_id"].n_unique() == compounds.height
    featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
    columns = [f"emb_{index:04d}" for index in range(mp.output_dim)]
    frames = {depth: [] for depth in depths}
    max_difference = 0.0
    with torch.inference_mode():
        for start in range(0, compounds.height, cfg["batch_size"]):
            batch = compounds.slice(start, cfg["batch_size"])
            graphs = [featurizer(Chem.MolFromSmiles(s)) for s in batch["standardized_smiles"]]
            graph = BatchMolGraph(graphs)
            graph.to(cfg["device"])
            depth_states = states(mp, graph, depths)
            reference = aggregation(mp(graph), graph.batch)
            for depth in depths:
                values = aggregation(depth_states[depth], graph.batch).cpu().numpy()
                assert values.shape == (batch.height, mp.output_dim)
                assert np.isfinite(values).all()
                frames[depth].append(
                    pl.DataFrame(values, schema=columns).insert_column(0, batch["compound_id"])
                )
                if depth == mp.depth:
                    max_difference = max(
                        max_difference, float(np.max(np.abs(values - reference.cpu().numpy())))
                    )
            if start % 1024 == 0:
                print(f"{start + batch.height}/{compounds.height}", flush=True)
    assert max_difference <= cfg["validation_atol"]
    for depth in depths:
        pl.concat(frames[depth]).write_parquet(out / f"depth-{depth}.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "depths": depths,
                "rows": compounds.height,
                "dimensions": mp.output_dim,
                "final_vs_chemprop_max_abs": max_difference,
                "definition": (
                    "finalize each directed-bond iteration with frozen W_o, "
                    "then mean-pool atoms"
                ),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
