"""Recover challenge embeddings from the completed selected PubChem training state."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop import nn as cnn
from chemprop.data import BatchMolGraph
from rdkit import Chem

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
ckpt = torch.load(next(inp.glob("*-best.pt")), map_location="cpu", weights_only=False)
c = pl.read_parquet(
    next(p for p in inp.iterdir() if p.name.partition("-")[2] == "compounds.parquet")
).sort("compound_id")
assert c.height == 6896 and c["compound_id"].n_unique() == 6896
torch.set_num_threads(4)
torch.backends.cuda.matmul.allow_tf32 = False
mp = cnn.BondMessagePassing(d_h=600, depth=5).cuda().eval()
head = (
    torch.nn.Sequential(
        torch.nn.Linear(600, 600), torch.nn.ReLU(), torch.nn.Dropout(0.1), torch.nn.Linear(600, 5)
    )
    .cuda()
    .eval()
)
for prefix, module in [("mp", mp), ("head", head)]:
    module.load_state_dict(
        {
            k.removeprefix(prefix + "."): v
            for k, v in ckpt["state_dict"].items()
            if k.startswith(prefix + ".")
        },
        strict=True,
    )
f = featurizers.SimpleMoleculeMolGraphFeaturizer()
values, preds = [], []
with torch.inference_mode():
    for start in range(0, c.height, 128):
        g = BatchMolGraph(
            [f(Chem.MolFromSmiles(s)) for s in c.slice(start, 128)["standardized_smiles"]]
        )
        g.to("cuda")
        emb = cnn.MeanAggregation()(mp(g), g.batch)
        values.append(emb.cpu().numpy())
        preds.append(head(emb).cpu().numpy() * ckpt["stds"] + ckpt["means"])
values = np.concatenate(values)
assert values.shape == (6896, 600) and np.isfinite(values).all()
pl.DataFrame(values, schema=[f"emb_{i:03d}" for i in range(600)]).insert_column(
    0, c["compound_id"]
).write_parquet(out / "pubchem600-embeddings.parquet")
pl.DataFrame(np.concatenate(preds), schema=ckpt["heads"]).insert_column(
    0, c["compound_id"]
).write_parquet(out / "challenge-pac50-predictions.parquet")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            rows=c.height,
            features=600,
            best_epoch=ckpt["epoch"],
            heads=ckpt["heads"],
            pretraining_config=ckpt["config"],
            source_run="20260911T153035269729Z-8d6bae93",
            recovery="Strict best-state load and exact challenge input; no retraining",
        ),
        indent=2,
    )
)
print(values.shape, ckpt["epoch"])
