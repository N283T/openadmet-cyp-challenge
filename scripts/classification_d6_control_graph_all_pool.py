"""Extend frozen TDI-only graph refits to all held-fold molecules."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop.data import BatchMolGraph
from finetune_external_multitask import MultitaskModel, encoder_architecture
from rdkit import Chem

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
assigned = (
    pl.read_parquet(inp / cfg["assignments"])
    .filter(pl.col("strategy") == cfg["split_id"])
    .select("compound_id", "fold", "group")
)
labels = (
    pl.read_parquet(inp / cfg["labels"])
    .filter(pl.col("endpoint") == cfg["endpoint"])
    .select("compound_id", pl.col("value").cast(pl.Int64).alias("y_true"))
)
pool = assigned.join(labels, on="compound_id", how="left", validate="1:1").sort("compound_id")
if pool.height != 6146 or pool["compound_id"].n_unique() != 6146:
    raise ValueError("Incomplete graph molecular pool")
ids = compounds["compound_id"].to_list()
lookup = {value: i for i, value in enumerate(ids)}
featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
graphs = [featurizer(Chem.MolFromSmiles(smiles)) for smiles in compounds["standardized_smiles"]]
checkpoint = torch.load(inp / cfg["encoder"], map_location="cpu", weights_only=False)
encoder = {
    key.removeprefix("mp."): value
    for key, value in checkpoint["state_dict"].items()
    if key.startswith("mp.")
}
architecture = encoder_architecture(checkpoint)
if "W_i.bias" not in encoder:
    architecture["bias"] = False
torch.set_num_threads(2)
parts = []
for fold, alias in enumerate(cfg["refits"]):
    held = pool.filter(pl.col("fold") == fold)
    positions = [lookup[value] for value in held["compound_id"]]
    model = MultitaskModel(encoder, 0.1, 2, architecture).to("cuda")
    state = torch.load(inp / alias, map_location="cpu", weights_only=False)
    model.load_state_dict(state["state_dict"])
    model.eval()
    predictions = []
    with torch.no_grad():
        for start in range(0, len(positions), 128):
            batch = BatchMolGraph([graphs[i] for i in positions[start : start + 128]])
            batch.to("cuda")
            predictions.append(model(batch).cpu().numpy()[:, 0])
    logits = np.concatenate(predictions)
    probability = 1 / (1 + np.exp(-logits))
    if not np.isfinite(probability).all():
        raise ValueError(f"Nonfinite graph probability in fold {fold}")
    parts.append(
        held.with_columns(
            pl.lit(cfg["endpoint"]).alias("endpoint"),
            pl.lit(cfg["split_id"]).alias("split_id"),
            pl.lit(cfg["model"]).alias("model"),
            pl.Series("score", probability),
        )
    )
    print(f"fold {fold} complete", flush=True)
pl.concat(parts).sort("compound_id").write_parquet(out / "all-pool-oof.parquet")
