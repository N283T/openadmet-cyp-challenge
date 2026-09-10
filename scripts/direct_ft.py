"""Direct CheMeleon fine-tuning, separate masked multitask regression and TDI."""

import gc
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop.data import BatchMolGraph
from rdkit import Chem
from split_baselines import metrics

from cyp_challenge.log2fc_pretrain import Network, partition, target_scaler
from cyp_challenge.neural_losses import masked_loss
from cyp_challenge.storage import digest, write_json


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    torch.set_num_threads(4)
    torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    device = torch.device("cuda")
    checkpoint_path = next(inp.glob("*-chemeleon_mp.pt"))
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    compounds = pl.read_parquet(next(inp.glob("*-compounds.parquet")))
    labels = pl.read_parquet(next(inp.glob("*-measurements.parquet")))
    assignments = pl.read_parquet(next(inp.glob("*-assignments.parquet"))).filter(
        pl.col("strategy") == "umap"
    )
    feat = featurizers.SimpleMoleculeMolGraphFeaturizer()
    graphs = {
        cid: feat(Chem.MolFromSmiles(smi))
        for cid, smi in compounds.select("compound_id", "standardized_smiles").iter_rows()
    }

    def graph(ids):
        g = BatchMolGraph([graphs[i] for i in ids])
        g.to(device)
        return g

    rows, scores = [], []
    for binary in [False, True]:
        torch.manual_seed(cfg["seed"] + cfg["fold"])
        track = "classification" if binary else "regression"
        heads = [ep for ep in cfg["endpoints"] if ep.endswith("is_TDI") == binary]
        wide = (
            labels.filter(pl.col("endpoint").is_in(heads))
            .pivot(on="endpoint", index="compound_id", values="value")
            .join(
                assignments.select("compound_id", "fold"),
                on="compound_id",
                how="left",
                validate="1:1",
            )
            .sort("compound_id")
        )
        assert wide["fold"].null_count() == 0
        ids = wide["compound_id"].to_list()
        y = wide.select(heads).to_numpy().astype(np.float32)
        train, val, held = partition(
            ids, wide["fold"].to_numpy(), cfg["fold"], cfg["seed"], cfg["inner_validation_fraction"]
        )
        means, stds = (
            (np.zeros(len(heads), dtype=np.float32), np.ones(len(heads), dtype=np.float32))
            if binary
            else target_scaler(y, train)
        )
        z = (y - means) / stds
        role = np.full(len(ids), "outer_holdout", dtype=object)
        role[train] = "train"
        role[val] = "inner_validation"
        pl.DataFrame(
            {"compound_id": ids, "role": role.tolist(), "fold": wide["fold"]}
        ).write_parquet(out / f"{track}-membership.parquet")
        meta = dict(
            config=cfg,
            heads=heads,
            means=means.tolist(),
            stds=stds.tolist(),
            checkpoint_sha256=digest(checkpoint_path),
            architecture=checkpoint["hyper_parameters"],
            train_rows=len(train),
            inner_validation_rows=len(val),
            outer_rows=len(held),
            loss="mean of per-head masked BCE"
            if binary
            else "mean of per-head masked standardized MSE",
        )
        write_json(out / f"{track}-metadata.json", meta)
        model = Network(checkpoint, cfg, len(heads)).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=cfg["learning_rate"])
        rng = np.random.default_rng(cfg["seed"])
        best = float("inf")
        stale = 0

        def predict(indices, model=model, ids=ids):
            model.eval()
            values = []
            with torch.inference_mode():
                for start in range(0, len(indices), cfg["batch_size"]):
                    batch = indices[start : start + cfg["batch_size"]]
                    values.append(model(graph([ids[i] for i in batch])).cpu())
            return torch.cat(values)

        for epoch in range(cfg["max_epochs"]):
            model.train()
            losses = []
            for batch in np.array_split(
                rng.permutation(train), int(np.ceil(len(train) / cfg["batch_size"]))
            ):
                optimizer.zero_grad(set_to_none=True)
                p = model(graph([ids[i] for i in batch]))
                loss = masked_loss(p, torch.tensor(z[batch], device=device), binary)
                assert torch.isfinite(loss)
                loss.backward()
                optimizer.step()
                losses.append(float(loss.detach()))
            value = float(masked_loss(predict(val), torch.tensor(z[val]), binary))
            assert np.isfinite(value)
            if value < best:
                best = value
                stale = 0
                torch.save(
                    dict(
                        state_dict=model.state_dict(),
                        metadata=meta,
                        epoch=epoch,
                        validation_loss=value,
                    ),
                    out / f"{track}-best.pt",
                )
            else:
                stale += 1
            torch.save(
                dict(
                    state_dict=model.state_dict(),
                    optimizer=optimizer.state_dict(),
                    epoch=epoch,
                    metadata=meta,
                ),
                out / f"{track}-last.pt",
            )
            entry = dict(
                track=track, epoch=epoch, train_loss=float(np.mean(losses)), val_loss=value
            )
            with (out / "history.jsonl").open("a") as f:
                f.write(json.dumps(entry) + "\n")
            print(json.dumps(entry), flush=True)
            if stale >= cfg["patience"]:
                break
        saved = torch.load(out / f"{track}-best.pt", map_location="cpu", weights_only=False)
        model.load_state_dict(saved["state_dict"])
        raw = predict(held)
        predictions = raw.sigmoid().numpy() if binary else raw.numpy() * stds + means
        for j, ep in enumerate(heads):
            mask = np.isfinite(y[held, j])
            selected = held[mask]
            p = predictions[mask, j]
            truth = (
                labels.filter(pl.col("endpoint") == ep)
                .join(
                    pl.DataFrame({"compound_id": [ids[i] for i in selected], "prediction": p}),
                    on="compound_id",
                    validate="1:1",
                )
                .sort("compound_id")
            )
            score = dict(
                endpoint=ep,
                fold=cfg["fold"],
                n=truth.height,
                **metrics(
                    truth["value"].to_numpy(),
                    truth["prediction"].to_numpy(),
                    truth["conf_low"].to_numpy(),
                    truth["conf_high"].to_numpy(),
                    binary,
                ),
            )
            scores.append(score)
            print(json.dumps(score), flush=True)
            rows.append(
                truth.select(
                    "compound_id", "prediction", pl.col("value").alias("y_true")
                ).with_columns(
                    pl.lit(ep).alias("endpoint"),
                    pl.lit(cfg["fold"]).alias("fold"),
                    pl.lit(cfg["split_run"]).alias("split_id"),
                    pl.lit("umap").alias("strategy"),
                )
            )
        del model, optimizer
        gc.collect()
        torch.cuda.empty_cache()
    pl.concat(rows).write_parquet(out / "oof.parquet")
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            metrics=scores,
            checkpoint_sha256=digest(checkpoint_path),
            torch=torch.__version__,
            device=torch.cuda.get_device_name(),
        ),
    )


if __name__ == "__main__":
    main()
