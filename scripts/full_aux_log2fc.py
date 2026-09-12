"""Refit adopted auxiliary on every measured compound and export cross-fitted train features."""

import json
import math
import os
import time
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop.data import BatchMolGraph
from rdkit import Chem

from cyp_challenge.log2fc_pretrain import Network, target_scaler
from cyp_challenge.storage import digest, write_json


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    reg = json.loads((inp / cfg["registry"]).read_text())
    recipe = reg["recipe"]
    summaries = [json.loads((inp / p).read_text()) for p in cfg["pair_summaries"]]
    epochs_list = [s["best_epoch"] + 1 for s in summaries]
    epochs = int(np.floor(np.median(epochs_list) + 0.5))
    assert len(epochs_list) == 10 and 1 <= epochs <= recipe["max_epochs"]
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    torch.manual_seed(recipe["seed"])
    torch.set_num_threads(3)
    torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
    nested = pl.read_parquet(inp / cfg["nested"])
    cf = (
        nested.filter(pl.col("fold") == pl.col("outer_fold"))
        .select("compound_id", "fold", *[f"aux_{i}" for i in range(4)])
        .sort("compound_id")
    )
    assert cf.height == cf["compound_id"].n_unique() == 6146
    targets = pl.read_parquet(inp / cfg["targets"])
    heads = sorted(targets["head"].unique())
    assert all(s["heads"] == heads for s in summaries)
    wide = (
        targets.pivot(on="head", index="compound_id", values="value")
        .join(
            compounds.select("compound_id", "standardized_smiles"), on="compound_id", validate="1:1"
        )
        .sort("compound_id")
    )
    assert wide.height == targets["compound_id"].n_unique() == 4376
    assert set(wide["compound_id"]) <= set(cf["compound_id"])
    y = wide.select(heads).to_numpy().astype(np.float32)
    means, stds = target_scaler(y, np.arange(wide.height))
    z = (y - means) / stds
    cp = inp / cfg["checkpoint"]
    assert digest(cp) == recipe["checkpoint_sha256"]
    checkpoint = torch.load(cp, map_location="cpu", weights_only=True)
    device = torch.device("cuda")
    model = Network(checkpoint, recipe, 4).to(device)
    optimizer = torch.optim.AdamW(
        [
            {"params": model.mp.parameters(), "lr_scale": recipe["encoder_lr_scale"]},
            {
                "params": list(model.bn.parameters()) + list(model.head.parameters()),
                "lr_scale": 1.0,
            },
        ],
        lr=recipe["learning_rate"],
        weight_decay=recipe["weight_decay"],
    )
    f = featurizers.SimpleMoleculeMolGraphFeaturizer()
    graphs = {
        cid: f(Chem.MolFromSmiles(s))
        for cid, s in compounds.select("compound_id", "standardized_smiles").iter_rows()
    }
    ids = wide["compound_id"].to_list()

    def graph(batch_ids):
        g = BatchMolGraph([graphs[c] for c in batch_ids])
        g.to(device)
        return g

    batch_size = recipe["batch_size"]
    steps = math.ceil(len(ids) / batch_size)
    warmup = recipe["warmup_epochs"] * steps
    total = recipe["max_epochs"] * steps
    peak = recipe["learning_rate"] * recipe["lr_ratio"]
    final_lr = recipe["learning_rate"] * 0.1
    rng, step = np.random.default_rng(recipe["seed"]), 0
    membership = wide.select("compound_id").with_columns(pl.lit("full_aux_train").alias("role"))
    membership.write_parquet(out / "membership.parquet")
    metadata = {
        "recipe": recipe,
        "heads": heads,
        "target_means": means.tolist(),
        "target_stds": stds.tolist(),
        "epochs": epochs,
        "cv_epoch_counts": epochs_list,
        "epoch_selection": "Rounded median of ten prior pair-model best epoch counts",
        "lr_horizon_epochs": recipe["max_epochs"],
        "train_compounds": len(ids),
        "train_labels": int(np.isfinite(y).sum()),
        "heldout_compounds": 0,
        "checkpoint_provenance": json.loads((inp / cfg["provenance"]).read_text()),
    }
    write_json(out / "metadata.json", metadata)
    for epoch in range(epochs):
        started = time.monotonic()
        model.train()
        loss_sum, count = 0.0, 0
        for batch in np.array_split(rng.permutation(len(ids)), steps):
            lr = (
                recipe["learning_rate"] + (peak - recipe["learning_rate"]) * step / max(1, warmup)
                if step < warmup
                else peak * (final_lr / peak) ** ((step - warmup) / max(1, total - warmup))
            )
            for group in optimizer.param_groups:
                group["lr"] = lr * group["lr_scale"]
            truth = torch.tensor(z[batch], device=device)
            mask = torch.isfinite(truth)
            optimizer.zero_grad(set_to_none=True)
            pred = model(graph([ids[i] for i in batch]))
            loss = ((pred[mask] - truth[mask]) ** 2).mean()
            assert torch.isfinite(loss)
            loss.backward()
            optimizer.step()
            loss_sum += float(loss.detach()) * int(mask.sum())
            count += int(mask.sum())
            step += 1
        row = {
            "epoch": epoch + 1,
            "train_mse_z": loss_sum / count,
            "learning_rate": lr,
            "seconds": time.monotonic() - started,
        }
        with (out / "history.jsonl").open("a") as stream:
            stream.write(json.dumps(row) + "\n")
        print(json.dumps(row), flush=True)
    torch.save(
        {
            "state_dict": model.cpu().state_dict(),
            "metadata": metadata,
            "architecture": checkpoint["hyper_parameters"],
            "epoch": epochs - 1,
        },
        out / "full.pt",
    )
    model.to(device).eval()
    test = compounds.join(cf.select("compound_id"), on="compound_id", how="anti").sort(
        "compound_id"
    )
    assert test.height == 750
    values = []
    with torch.inference_mode():
        for start in range(0, test.height, batch_size):
            batch = test.slice(start, batch_size)
            values.append(model(graph(batch["compound_id"].to_list())).cpu().numpy() * stds + means)
    p = np.concatenate(values)
    assert np.isfinite(p).all()
    test_aux = test.select("compound_id").with_columns(
        pl.lit(-1, dtype=cf.schema["fold"]).alias("fold"),
        *[pl.Series(f"aux_{i}", p[:, i]) for i in range(4)],
    )
    result = pl.concat(
        [
            cf.with_columns(pl.lit("cross_fitted_train").alias("role")),
            test_aux.with_columns(pl.lit("full_fit_test").alias("role")),
        ]
    ).sort("compound_id")
    assert result.height == result["compound_id"].n_unique() == 6896
    result.write_parquet(out / "full-aux.parquet")
    write_json(
        out / "summary.json",
        {
            **metadata,
            "feature_rows": result.height,
            "test_rows": test.height,
            "feature_heads": dict(zip([f"aux_{i}" for i in range(4)], heads, strict=True)),
            "train_features": "Four-model mean excluding recipient fold, including anchors",
            "full_checkpoint_sha256": digest(out / "full.pt"),
            "feature_sha256": digest(out / "full-aux.parquet"),
        },
    )


if __name__ == "__main__":
    main()
