"""External-only Chemprop initialization with honest classification UMAP inner/refit."""

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
from sklearn.metrics import matthews_corrcoef, roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


compounds, labels, assignments, testmap, aux = [pl.read_parquet(path(i)) for i in range(5)]
assignments = assignments.filter(pl.col("strategy") == "umap")
assert assignments.height == 6146
class_heads = ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]
class_labels = labels.filter(pl.col("endpoint").is_in(class_heads))
assert len(class_labels) == 5081
aux = aux.filter(
    (pl.col("kind") == "logfc")
    | ((pl.col("kind") == "tdi") & pl.col("enzyme").is_in(["CYP2D6", "CYP3A4"]))
).select("compound_id", (pl.col("enzyme") + "_" + pl.col("kind")).alias("endpoint"), "value")
direct = labels.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")).select(
    "compound_id", "endpoint", "value"
)
aux = pl.concat([aux, direct])
aux_heads = sorted(aux["endpoint"].unique())
assert len(aux_heads) == 10
heads = class_heads + aux_heads
values = pl.concat([class_labels.select("compound_id", "endpoint", "value"), aux])
assert not values.select("compound_id", "endpoint").is_duplicated().any()
table = (
    compounds.sort("compound_id")
    .join(
        assignments.select("compound_id", "fold"),
        on="compound_id",
        how="left",
        validate="1:1",
        maintain_order="left",
    )
    .join(
        values.pivot(on="endpoint", index="compound_id", values="value"),
        on="compound_id",
        how="left",
        validate="1:1",
        maintain_order="left",
    )
)
assert len(table) == 6896
ids = table["compound_id"].to_numpy()
y = table.select(heads).to_numpy().astype(np.float32)
folds = table["fold"].fill_null(-1).to_numpy()
test_idx = np.flatnonzero(np.isin(ids, testmap["compound_id"].to_numpy()))
assert len(test_idx) == 750 and np.all(folds[test_idx] == -1)
checkpoint = torch.load(path(5), map_location="cpu", weights_only=False)
state = {
    k.removeprefix("mp."): v for k, v in checkpoint["state_dict"].items() if k.startswith("mp.")
}
architecture = encoder_architecture(checkpoint)
torch.set_num_threads(2)
device = "cuda"
f = featurizers.SimpleMoleculeMolGraphFeaturizer()
graphs = [f(Chem.MolFromSmiles(s)) for s in table["standardized_smiles"]]


def graph(ix):
    b = BatchMolGraph([graphs[int(i)] for i in ix])
    b.to(device)
    return b


def predict(model, ix):
    model.eval()
    with torch.no_grad():
        return np.concatenate(
            [model(graph(ix[j : j + 128])).cpu().numpy() for j in range(0, len(ix), 128)]
        )


def bce(pred, truth):
    mask = torch.isfinite(truth)
    err = torch.nn.functional.binary_cross_entropy_with_logits(
        pred, torch.nan_to_num(truth), reduction="none"
    )
    present = mask.any(0)
    return ((err * mask).sum(0) / mask.sum(0).clamp_min(1))[present].mean()


all_oof, all_inner, all_test, summaries = [], [], [], []
(out / "target-contract.json").write_text(
    json.dumps(
        dict(
            heads=heads,
            counts=[int(np.isfinite(y[:, j]).sum()) for j in range(12)],
            architecture=architecture,
            config=cfg,
        ),
        indent=2,
    )
)
for recipe in ["tdi_only", "tdi_aux"]:
    nheads = 2 if recipe == "tdi_only" else 12
    relevant = np.isfinite(y[:, :nheads]).any(1)
    for outer in range(5):
        inner = (outer + 1) % 5
        directory = out / f"{recipe}-fold{outer}"
        directory.mkdir()
        outer_ix = np.flatnonzero(folds == outer)
        inner_ix = np.flatnonzero((folds == inner) & np.isfinite(y[:, :2]).any(1))
        selected = None
        best_state = None
        thresholds = None
        for phase in ["inner", "refit"]:
            trainmask = (folds >= 0) & (folds != outer) & relevant
            if phase == "inner":
                trainmask &= folds != inner
            train_ix = np.flatnonzero(trainmask)
            assert not set(train_ix) & set(outer_ix)
            if phase == "inner":
                assert not set(train_ix) & set(inner_ix)
            mean = np.zeros(nheads, dtype=np.float32)
            std = np.ones(nheads, dtype=np.float32)
            if nheads > 2:
                mean[2:] = np.nanmean(y[train_ix, 2:nheads], axis=0)
                std[2:] = np.maximum(np.nanstd(y[train_ix, 2:nheads], axis=0), 1e-6)
            z = (y[:, :nheads] - mean) / std
            assert np.isfinite(mean).all() and np.isfinite(std).all()
            (directory / f"{phase}-scaling.json").write_text(
                json.dumps(
                    dict(
                        mean=mean.tolist(),
                        std=std.tolist(),
                        train_counts=np.isfinite(z[train_ix]).sum(0).tolist(),
                    )
                )
            )
            table.select("compound_id", "fold").with_columns(
                pl.Series("training", trainmask),
                pl.Series("outer_holdout", folds == outer),
                pl.Series("inner_validation", (folds == inner) & (phase == "inner")),
                pl.lit(outer).alias("outer_fold"),
                pl.lit(phase).alias("phase"),
            ).write_parquet(directory / f"{phase}-membership.parquet")
            torch.manual_seed(42 + outer)
            np_rng = np.random.default_rng(42 + outer)
            model = MultitaskModel(state, 0.1, nheads, architecture).to(device)
            optimizer = torch.optim.AdamW(
                [
                    {"params": model.mp.parameters(), "lr": 1e-5},
                    {"params": model.head.parameters(), "lr": 1e-4},
                ],
                weight_decay=1e-4,
            )
            best, stale = float("inf"), 0
            limit = 80 if phase == "inner" else selected
            logs = []
            for epoch in range(1, limit + 1):
                model.train()
                losses = []
                for batch in np.array_split(
                    np_rng.permutation(train_ix), int(np.ceil(len(train_ix) / 64))
                ):
                    truth = torch.tensor(z[batch], device=device)
                    pred = model(graph(batch))
                    loss = (
                        bce(pred[:, :2], truth[:, :2])
                        if torch.isfinite(truth[:, :2]).any()
                        else pred[:, :2].sum() * 0
                    )
                    if nheads > 2:
                        mask = torch.isfinite(truth[:, 2:])
                        err = (pred[:, 2:] - torch.nan_to_num(truth[:, 2:])).square()
                        present = mask.any(0)
                        if present.any():
                            loss = (
                                loss
                                + 0.1
                                * ((err * mask).sum(0) / mask.sum(0).clamp_min(1))[present].mean()
                            )
                    optimizer.zero_grad()
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(model.parameters(), 5)
                    optimizer.step()
                    losses.append(float(loss.detach()))
                metric = None
                if phase == "inner":
                    p = predict(model, inner_ix)[:, :2]
                    metric = float(bce(torch.tensor(p), torch.tensor(y[inner_ix, :2])))
                    if metric < best - 1e-6:
                        best, selected, stale = metric, epoch, 0
                        best_state = copy.deepcopy(
                            {k: v.cpu() for k, v in model.state_dict().items()}
                        )
                    else:
                        stale += 1
                logs.append(dict(epoch=epoch, loss=float(np.mean(losses)), inner_BCE=metric))
                (directory / f"{phase}-history.json").write_text(json.dumps(logs))
                if epoch % 10 == 0:
                    print(recipe, outer, phase, epoch, metric, flush=True)
                if phase == "inner" and stale >= 10:
                    break
            if phase == "inner":
                model.load_state_dict(best_state)
                inner_p = 1 / (1 + np.exp(-predict(model, inner_ix)[:, :2]))
                thresholds = []
                for j, ep in enumerate(class_heads):
                    mask = np.isfinite(y[inner_ix, j])
                    yy = y[inner_ix, j][mask]
                    pp = inner_p[mask, j]
                    candidates = np.unique(np.r_[0.5, np.linspace(0.05, 0.95, 91)])
                    scores = [matthews_corrcoef(yy, pp >= t) for t in candidates]
                    threshold = float(candidates[np.argmax(scores)])
                    thresholds.append(threshold)
                    all_inner.append(
                        pl.DataFrame(
                            dict(
                                compound_id=ids[inner_ix][mask],
                                endpoint=ep,
                                fold=inner,
                                outer_fold=outer,
                                y_true=yy,
                                probability=pp,
                                threshold=threshold,
                                method=recipe,
                            )
                        )
                    )
                torch.save(
                    dict(
                        state_dict=best_state,
                        epoch=selected,
                        heads=heads[:nheads],
                        thresholds=thresholds,
                    ),
                    directory / "inner-best.pt",
                )
            else:
                torch.save(
                    dict(
                        state_dict=model.cpu().state_dict(),
                        epochs=selected,
                        heads=heads[:nheads],
                        architecture=architecture,
                    ),
                    directory / "refit.pt",
                )
                model.to(device)
                for ix, collection, is_test in [
                    (outer_ix, all_oof, False),
                    (test_idx, all_test, True),
                ]:
                    pp = 1 / (1 + np.exp(-predict(model, ix)[:, :2]))
                    for j, ep in enumerate(class_heads):
                        keep = np.ones(len(ix), bool) if is_test else np.isfinite(y[ix, j])
                        collection.append(
                            pl.DataFrame(
                                dict(
                                    compound_id=ids[ix][keep],
                                    endpoint=ep,
                                    fold=outer,
                                    split_id="20260909T135130268493Z-edba6fc0",
                                    y_true=y[ix, j][keep],
                                    probability=pp[keep, j],
                                    prediction=(pp[keep, j] >= thresholds[j]).astype(int),
                                    threshold=thresholds[j],
                                    method=recipe,
                                )
                            )
                        )
        summaries.append(
            dict(method=recipe, fold=outer, selected_epochs=selected, thresholds=thresholds)
        )
        pl.concat(all_oof).write_parquet(out / "classification-oof.parquet")
        pl.concat(all_inner).write_parquet(out / "inner-probabilities.parquet")
        pl.concat(all_test).write_parquet(out / "fold-test-probabilities.parquet")
        (out / "progress.json").write_text(json.dumps(summaries, indent=2))
        print("DONE", recipe, outer, selected, flush=True)
result = pl.concat(all_oof)
metrics = []
for (method, ep), t in result.group_by("method", "endpoint"):
    metrics.append(
        dict(
            method=method,
            endpoint=ep,
            n=t.height,
            MCC=float(matthews_corrcoef(t["y_true"], t["prediction"])),
            ROC_AUC=float(roc_auc_score(t["y_true"], t["probability"])),
        )
    )
pl.DataFrame(metrics).write_csv(out / "metrics.csv")
pl.concat(all_test).group_by("compound_id", "endpoint", "method").agg(
    pl.col("probability").mean(), pl.col("threshold").mean()
).with_columns(
    (pl.col("probability") >= pl.col("threshold")).cast(pl.Int8).alias("prediction")
).write_parquet(out / "test-probabilities.parquet")
print(pl.DataFrame(metrics).write_csv())
