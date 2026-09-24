"""Refit the TDI-only graph control on scaffold folds with inner epoch selection."""

import copy
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
if assigned.height != 6146 or assigned["compound_id"].n_unique() != 6146:
    raise ValueError("Incomplete scaffold molecular pool")
labels = pl.read_parquet(inp / cfg["labels"]).filter(
    pl.col("endpoint").is_in(["CYP2D6_is_TDI", "CYP3A4_is_TDI"])
)
wide = labels.select("compound_id", "endpoint", "value").pivot(
    on="endpoint", index="compound_id", values="value"
)
table = compounds.join(assigned, on="compound_id", how="left", validate="1:1").join(
    wide, on="compound_id", how="left", validate="1:1"
)
folds = table["fold"].fill_null(-1).to_numpy()
y = table.select("CYP2D6_is_TDI", "CYP3A4_is_TDI").to_numpy().astype(np.float32)
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
featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
graphs = [featurizer(Chem.MolFromSmiles(smiles)) for smiles in table["standardized_smiles"]]


def graph(indices: np.ndarray) -> BatchMolGraph:
    batch = BatchMolGraph([graphs[int(i)] for i in indices])
    batch.to("cuda")
    return batch


def predict(model: MultitaskModel, indices: np.ndarray) -> np.ndarray:
    model.eval()
    with torch.no_grad():
        return np.concatenate(
            [
                model(graph(indices[start : start + 128])).cpu().numpy()
                for start in range(0, len(indices), 128)
            ]
        )


def bce(logits: torch.Tensor, truth: torch.Tensor) -> torch.Tensor:
    mask = torch.isfinite(truth)
    error = torch.nn.functional.binary_cross_entropy_with_logits(
        logits, torch.nan_to_num(truth), reduction="none"
    )
    present = mask.any(0)
    return ((error * mask).sum(0) / mask.sum(0).clamp_min(1))[present].mean()


parts = []
epochs = []
relevant = np.isfinite(y).any(1)
for outer in range(5):
    inner = (outer + 1) % 5
    held = np.flatnonzero(folds == outer)
    inner_validation = np.flatnonzero((folds == inner) & relevant)
    best_epoch = None
    for phase in ("inner", "refit"):
        trainmask = (folds >= 0) & (folds != outer) & relevant
        if phase == "inner":
            trainmask &= folds != inner
        train_indices = np.flatnonzero(trainmask)
        torch.manual_seed(42 + outer)
        rng = np.random.default_rng(42 + outer)
        model = MultitaskModel(encoder, 0.1, 2, architecture).to("cuda")
        optimizer = torch.optim.AdamW(
            [
                {"params": model.mp.parameters(), "lr": 1e-5},
                {"params": model.head.parameters(), "lr": 1e-4},
            ],
            weight_decay=1e-4,
        )
        best, stale = float("inf"), 0
        limit = 80 if phase == "inner" else best_epoch
        for epoch in range(1, limit + 1):
            model.train()
            for batch in np.array_split(
                rng.permutation(train_indices), int(np.ceil(len(train_indices) / 64))
            ):
                truth = torch.tensor(y[batch], device="cuda")
                loss = bce(model(graph(batch)), truth)
                optimizer.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 5)
                optimizer.step()
            if phase == "inner":
                metric = float(
                    bce(
                        torch.tensor(predict(model, inner_validation)),
                        torch.tensor(y[inner_validation]),
                    )
                )
                if metric < best - 1e-6:
                    best, best_epoch, stale = metric, epoch, 0
                    best_state = copy.deepcopy(
                        {key: value.cpu() for key, value in model.state_dict().items()}
                    )
                else:
                    stale += 1
                if stale >= 10:
                    break
        if phase == "inner":
            model.load_state_dict(best_state)
        else:
            logits = predict(model, held)[:, 0]
            probability = 1 / (1 + np.exp(-logits))
            if not np.isfinite(probability).all():
                raise ValueError(f"Nonfinite scaffold graph probability: {outer}")
            held_ids = table["compound_id"].to_numpy()[held]
            lookup = {value: index for index, value in enumerate(held_ids)}
            pool = (
                assigned.filter(pl.col("fold") == outer)
                .join(
                    labels.filter(pl.col("endpoint") == cfg["endpoint"]).select(
                        "compound_id", pl.col("value").cast(pl.Int64).alias("y_true")
                    ),
                    on="compound_id",
                    how="left",
                    validate="1:1",
                )
                .sort("compound_id")
            )
            order = [lookup[value] for value in pool["compound_id"]]
            parts.append(
                pool.with_columns(
                    pl.lit(cfg["endpoint"]).alias("endpoint"),
                    pl.lit(cfg["split_id"]).alias("split_id"),
                    pl.lit(cfg["model"]).alias("model"),
                    pl.Series("score", probability[order]),
                )
            )
            epochs.append({"fold": outer, "selected_epochs": best_epoch})
    print(f"fold {outer} epochs {best_epoch} complete", flush=True)
pl.concat(parts).sort("compound_id").write_parquet(out / "all-pool-oof.parquet")
(out / "selected-epochs.json").write_text(json.dumps(epochs, indent=2) + "\n")
