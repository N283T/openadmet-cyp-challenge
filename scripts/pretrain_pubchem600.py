"""Scratch 600-wide five-head PubChem pAC50 encoder; no challenge labels."""

import copy
import json
import math
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop import nn as cnn
from chemprop.data import BatchMolGraph
from rdkit import Chem
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(n):
    return pl.read_parquet(
        next(p for p in inp.iterdir() if p.name.partition("-")[2] == n + ".parquet")
    )


c = read("external-compounds").sort("compound_id")
t = read("targets")
heads = sorted(t["isoform"].unique())
h = c.join(
    t.pivot(on="isoform", index="compound_id", values="value"), on="compound_id", validate="1:1"
).sort("compound_id")
y = h.select(heads).to_numpy().astype(np.float32)
rng = np.random.default_rng(42)
keys = rng.permutation(sorted(h["connectivity"].unique()))
val_keys = set(keys[: max(1, len(keys) // 10)])
is_val = np.array([k in val_keys for k in h["connectivity"]])
train, validation = np.flatnonzero(~is_val), np.flatnonzero(is_val)
means, stds = np.nanmean(y[train], axis=0), np.nanstd(y[train], axis=0)
z = (y - means) / stds
h.select("compound_id", "connectivity").with_columns(
    pl.Series("role", np.where(is_val, "validation", "training"))
).write_parquet(out / "membership.parquet")
torch.manual_seed(42)
torch.set_num_threads(4)
torch.backends.cuda.matmul.allow_tf32 = False
featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
graphs = [featurizer(Chem.MolFromSmiles(s)) for s in h["standardized_smiles"]]


def graph(indices):
    g = BatchMolGraph([graphs[int(i)] for i in indices])
    g.to("cuda")
    return g


class Model(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.mp = cnn.BondMessagePassing(d_h=600, depth=5)
        self.head = torch.nn.Sequential(
            torch.nn.Linear(600, 600),
            torch.nn.ReLU(),
            torch.nn.Dropout(0.1),
            torch.nn.Linear(600, 5),
        )

    def embed(self, g):
        return cnn.MeanAggregation()(self.mp(g), g.batch)

    def forward(self, g):
        return self.head(self.embed(g))


model = Model().cuda()
opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
best, stale, best_epoch = float("inf"), 0, -1
history = []
for epoch in range(cfg["max_epochs"]):
    lr = (
        1e-3 * (epoch + 1) / 3
        if epoch < 3
        else 1e-4
        + 0.5
        * (1e-3 - 1e-4)
        * (1 + math.cos(math.pi * (epoch - 3) / max(1, cfg["max_epochs"] - 4)))
    )
    for pg in opt.param_groups:
        pg["lr"] = lr
    model.train()
    total, count = 0.0, 0
    for batch in np.array_split(rng.permutation(train), int(np.ceil(len(train) / 128))):
        truth = torch.as_tensor(z[batch], device="cuda")
        mask = torch.isfinite(truth)
        opt.zero_grad(set_to_none=True)
        pred = model(graph(batch))
        squared = (pred - torch.nan_to_num(truth)) ** 2
        per_head = (squared * mask).sum(0) / mask.sum(0).clamp_min(1)
        loss = per_head[mask.any(0)].mean()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5)
        opt.step()
        total += float(loss.detach()) * len(batch)
        count += len(batch)
    model.eval()
    with torch.inference_mode():
        pv = np.concatenate(
            [
                model(graph(batch)).cpu().numpy()
                for batch in np.array_split(validation, int(np.ceil(len(validation) / 128)))
            ]
        )
    vmse = float(np.nanmean(np.nanmean((pv - z[validation]) ** 2, axis=0)))
    assert np.isfinite(vmse)
    row = dict(epoch=epoch, train_mse=total / count, validation_mse=vmse, lr=lr)
    history.append(row)
    if vmse < best:
        best, best_epoch, stale = vmse, epoch, 0
        state = copy.deepcopy(model.state_dict())
        torch.save(
            dict(
                state_dict=state,
                optimizer=opt.state_dict(),
                epoch=epoch,
                heads=heads,
                means=means,
                stds=stds,
                config=cfg,
            ),
            out / "best.pt",
        )
    else:
        stale += 1
    torch.save(
        dict(
            state_dict=model.state_dict(),
            optimizer=opt.state_dict(),
            epoch=epoch,
            heads=heads,
            means=means,
            stds=stds,
        ),
        out / "last.pt",
    )
    with (out / "history.jsonl").open("a") as stream:
        stream.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)
    if stale >= cfg["patience"]:
        break
model.load_state_dict(state, strict=True)
model.eval()
with torch.inference_mode():
    pv = (
        np.concatenate(
            [
                model(graph(batch)).cpu().numpy()
                for batch in np.array_split(validation, int(np.ceil(len(validation) / 128)))
            ]
        )
        * stds
        + means
    )
metrics = []
for j, head in enumerate(heads):
    mask = np.isfinite(y[validation, j])
    yy, pp = y[validation, j][mask], pv[:, j][mask]
    metrics.append(
        dict(
            head=head,
            n=int(mask.sum()),
            rmse=float(np.sqrt(np.mean((yy - pp) ** 2))),
            spearman=float(spearmanr(yy, pp).statistic),
        )
    )
pl.DataFrame(metrics).write_csv(out / "external-validation-metrics.csv")
challenge = read("compounds").sort("compound_id")
values, ps = [], []
with torch.inference_mode():
    for start in range(0, challenge.height, 128):
        g = BatchMolGraph(
            [
                featurizer(Chem.MolFromSmiles(s))
                for s in challenge.slice(start, 128)["standardized_smiles"]
            ]
        )
        g.to("cuda")
        emb = model.embed(g)
        values.append(emb.cpu().numpy())
        ps.append(model.head(emb).cpu().numpy() * stds + means)
values = np.concatenate(values)
assert values.shape == (6896, 600) and np.isfinite(values).all()
pl.DataFrame(values, schema=[f"emb_{i:03d}" for i in range(600)]).insert_column(
    0, challenge["compound_id"]
).write_parquet(out / "pubchem600-embeddings.parquet")
pl.DataFrame(np.concatenate(ps), schema=heads).insert_column(
    0, challenge["compound_id"]
).write_parquet(out / "challenge-pac50-predictions.parquet")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            training_compounds=len(train),
            validation_compounds=len(validation),
            best_epoch=best_epoch,
            epochs=epoch + 1,
            best_validation_mse=best,
            heads=heads,
            validation_metrics=metrics,
            challenge_labels_used=False,
        ),
        indent=2,
    )
)
