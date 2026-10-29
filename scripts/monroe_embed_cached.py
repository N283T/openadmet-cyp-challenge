"""Frozen MONROE embeddings from cached conformer graphs and an arbitrary checkpoint.

Reads the graph cache written by monroe_graph_cache.py, so no conformers are generated.
The embedding is the mean over all cached conformers, matching the fine-tuning readout.
"""

import json
import os
import pickle
from pathlib import Path

import numpy as np
import polars as pl
import torch
from monroe.eval.embed import to_pyg
from monroe.model.ckpt import load_ckpt
from torch_geometric.data import Batch


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    device = cfg["device"]
    checkpoint = Path(os.environ["TMPDIR"]) / "checkpoint"
    checkpoint.mkdir(parents=True, exist_ok=True)
    for name, alias in cfg["checkpoint"].items():
        if not (checkpoint / name).exists():
            (checkpoint / name).symlink_to(inp / alias)
    graphs = []
    for name in cfg["graph_files"]:
        with open(inp / name, "rb") as handle:
            graphs.append(pickle.load(handle))
    ids = sorted(graphs[0])
    encoder = load_ckpt(str(checkpoint)).to(device).eval()
    total = np.zeros((len(ids), cfg["embedding_dim"]), dtype=np.float64)
    with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16):
        for table in graphs:
            for start in range(0, len(ids), cfg["batch_size"]):
                chunk = ids[start : start + cfg["batch_size"]]
                batch = Batch.from_data_list([to_pyg(table[c]) for c in chunk]).to(device)
                total[start : start + len(chunk)] += encoder(batch)[0].float().cpu().numpy()
    total /= len(graphs)
    assert np.isfinite(total).all()
    columns = {f"emb_{i:03d}": total[:, i] for i in range(total.shape[1])}
    pl.DataFrame({"compound_id": ids, **columns}).write_parquet(out / "external-embeddings.parquet")
    (out / "summary.json").write_text(json.dumps({"compounds": len(ids), "draws": len(graphs)}))


if __name__ == "__main__":
    main()
