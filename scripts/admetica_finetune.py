"""Matched-CYP pIC50 fine-tuning and downstream TabPFN on adapted embeddings."""

import copy
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
from tabpfn import TabPFNRegressor
from tabpfn.constants import ModelVersion

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(n):
    return pl.read_parquet(next(inp.glob("*-" + n + ".parquet")))


labels = read("measurements")
a = read("parent-assignments")
c = read("compounds")
torch.set_num_threads(4)
torch.backends.cuda.matmul.allow_tf32 = False
featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
rows, audits = [], []


class Model(torch.nn.Module):
    def __init__(self, checkpoint):
        super().__init__()
        self.mp = cnn.BondMessagePassing(**checkpoint["hyper_parameters"])
        self.mp.load_state_dict(checkpoint["state_dict"], strict=True)
        self.head = torch.nn.Sequential(
            torch.nn.Linear(300, 256),
            torch.nn.ReLU(),
            torch.nn.Dropout(0.1),
            torch.nn.Linear(256, 1),
        )

    def embed(self, g):
        return cnn.MeanAggregation()(self.mp(g), g.batch)

    def forward(self, g):
        return self.head(self.embed(g)).squeeze(-1)


for path in sorted(inp.glob("*-encoder.pt")):
    enzyme = next(e for e in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"] if e.lower() in path.name)
    ep = enzyme + "_pIC50_direct_inhibition"
    with torch.serialization.safe_globals([cnn.Activation]):
        checkpoint = torch.load(path, map_location="cpu", weights_only=True)
    h = (
        labels.filter(pl.col("endpoint") == ep)
        .join(a, on="compound_id", validate="1:1")
        .join(c, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    graphs = [featurizer(Chem.MolFromSmiles(s)) for s in h["standardized_smiles"]]
    y = h["value"].to_numpy().astype(np.float32)

    def graph(indices, graphs=graphs):
        g = BatchMolGraph([graphs[int(i)] for i in indices])
        g.to("cuda")
        return g

    for fold in range(5):
        seed = 42 + fold
        torch.manual_seed(seed)
        rng = np.random.default_rng(seed)
        outer = (h["fold"].to_numpy() == fold) & ~h["is_anchor"].to_numpy()
        eligible = np.flatnonzero(~outer)
        order = rng.permutation(eligible)
        n_val = max(1, int(len(order) * 0.1))
        validation, training = order[:n_val], order[n_val:]
        mean, std = float(y[training].mean()), float(y[training].std())
        z = (y - mean) / max(std, 1e-6)
        model = Model(checkpoint).cuda()
        opt = torch.optim.AdamW(
            [
                {"params": model.mp.parameters(), "lr": cfg["encoder_lr"]},
                {"params": model.head.parameters(), "lr": cfg["head_lr"]},
            ],
            weight_decay=cfg["weight_decay"],
        )
        best, best_epoch, stale = float("inf"), -1, 0
        history = []
        directory = out / f"{enzyme}-fold{fold}"
        directory.mkdir()
        h.select("compound_id").with_columns(
            pl.Series(
                "role",
                [
                    "outer"
                    if outer[i]
                    else "inner_validation"
                    if i in set(validation)
                    else "training"
                    for i in range(h.height)
                ],
            ),
            pl.lit(ep).alias("endpoint"),
            pl.lit(fold).alias("fold"),
        ).write_parquet(directory / "membership.parquet")
        for epoch in range(cfg["max_epochs"]):
            model.train()
            losses = []
            for batch in np.array_split(
                rng.permutation(training), int(np.ceil(len(training) / 64))
            ):
                opt.zero_grad(set_to_none=True)
                loss = torch.nn.functional.mse_loss(
                    model(graph(batch)), torch.as_tensor(z[batch], device="cuda")
                )
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 5)
                opt.step()
                losses.append((len(batch), float(loss.detach())))
            model.eval()
            with torch.inference_mode():
                pv = np.concatenate(
                    [
                        model(graph(batch)).cpu().numpy()
                        for batch in np.array_split(validation, int(np.ceil(len(validation) / 128)))
                    ]
                )
            val = float(np.mean((pv - z[validation]) ** 2))
            assert np.isfinite(val)
            history.append(
                dict(
                    epoch=epoch,
                    train_mse=sum(n * loss_value for n, loss_value in losses) / len(training),
                    validation_mse=val,
                )
            )
            if val < best:
                best, best_epoch, stale = val, epoch, 0
                state = copy.deepcopy(model.state_dict())
                torch.save(
                    dict(
                        state_dict=state,
                        optimizer=opt.state_dict(),
                        epoch=epoch,
                        mean=mean,
                        std=std,
                        encoder_parameters=checkpoint["hyper_parameters"],
                        config=cfg,
                    ),
                    directory / "best.pt",
                )
            else:
                stale += 1
            if stale >= cfg["patience"]:
                break
        torch.save(
            dict(
                state_dict=model.state_dict(),
                optimizer=opt.state_dict(),
                epoch=epoch,
                mean=mean,
                std=std,
            ),
            directory / "last.pt",
        )
        pl.DataFrame(history).write_csv(directory / "history.csv")
        model.load_state_dict(state, strict=True)
        model.eval()
        with torch.inference_mode():
            emb = np.concatenate(
                [
                    model.embed(graph(batch)).cpu().numpy()
                    for batch in np.array_split(np.arange(h.height), int(np.ceil(h.height / 128)))
                ]
            )
            direct = (
                model.head(torch.as_tensor(emb[outer], device="cuda")).squeeze(-1).cpu().numpy()
                * std
                + mean
            )
        pl.DataFrame(emb, schema=[f"emb_{i:03d}" for i in range(300)]).insert_column(
            0, h["compound_id"]
        ).write_parquet(directory / "embeddings.parquet")
        tab = TabPFNRegressor.create_default_for_version(
            ModelVersion.V2_6,
            n_estimators=8,
            model_path=str(inp / cfg["weights"]),
            device="cuda",
            random_state=seed,
            auto_scale_n_estimators=False,
            n_preprocessing_jobs=1,
            memory_saving_mode="auto",
            show_progress_bar=False,
            ignore_pretraining_limits=True,
            fit_mode="fit_preprocessors",
        )
        tab.fit(emb[~outer], y[~outer])
        tp = tab.predict(emb[outer])
        tab.save_fit_state(directory / "tabpfn.tabpfn_fit")
        for name, pred in [("admetica_ft_direct", direct), ("admetica_ft_emb_tabpfn", tp)]:
            assert np.isfinite(pred).all()
            rows.append(
                h.filter(pl.Series(outer))
                .select(
                    "compound_id",
                    "endpoint",
                    "fold",
                    "is_neighbor",
                    "conf_low",
                    "conf_high",
                    pl.col("value").alias("y_true"),
                )
                .with_columns(
                    pl.lit(name).alias("method"),
                    pl.lit(cfg["split_run"]).alias("split_id"),
                    pl.Series("prediction", pred.astype(float)),
                )
            )
        audits.append(
            dict(
                endpoint=ep,
                fold=fold,
                training=len(training),
                inner_validation=len(validation),
                outer=int(outer.sum()),
                best_epoch=best_epoch,
                epochs=epoch + 1,
                best_validation_mse=best,
            )
        )
        print(json.dumps(audits[-1]), flush=True)
pl.concat(rows).write_parquet(out / "finetune-oof.parquet")
(out / "summary.json").write_text(json.dumps(dict(config=cfg, fits=audits), indent=2))
