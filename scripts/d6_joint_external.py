"""Joint D-MPNN fine-tuning on D6 potency and separate external efficacy heads."""

import copy
import importlib.util
import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop.data import BatchMolGraph
from d6_specialist_fit import MultitaskModel, encoder_architecture
from rdkit import Chem

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


def read(i):
    return pl.read_parquet(path(i))


refit_epochs = (
    {(r["fold"], r["method"]): r["best_epoch"] for r in json.loads(path(7).read_text())}
    if cfg.get("outer_refit")
    else None
)


compounds, labels, assignments, external, nested = [read(i) for i in range(5)]
compounds = compounds.sort("compound_id")
external = external.sort("compound_id")
assert compounds.height == 6896
forbidden = {
    Chem.MolToInchiKey(Chem.MolFromSmiles(s)).split("-")[0]
    for s in compounds["standardized_smiles"]
}
assert set(external["connectivity"]).isdisjoint(forbidden)
assert set(compounds["compound_id"]).isdisjoint(set(external["compound_id"]))
labels = labels.filter(pl.col("endpoint") == "CYP2D6_pIC50_direct_inhibition")
spec = importlib.util.spec_from_file_location("official", path(6))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
checkpoint = torch.load(path(5), map_location="cpu", weights_only=False)
state = {
    k.removeprefix("mp."): v for k, v in checkpoint["state_dict"].items() if k.startswith("mp.")
}
arch = encoder_architecture(checkpoint)
torch.set_num_threads(2)
torch.backends.cuda.matmul.allow_tf32 = False
device = "cuda"
allc = pl.concat(
    [
        compounds.select("compound_id", "standardized_smiles"),
        external.select("compound_id", "standardized_smiles"),
    ]
)
lookup = {c: i for i, c in enumerate(allc["compound_id"])}
n = allc.height
truth = np.full((n, 4), np.nan, dtype=np.float32)
for r in labels.iter_rows(named=True):
    truth[lookup[r["compound_id"]], 0] = r["value"]
truth[6896:, 1:] = (
    external.select("Max_Response", "Activity at 57.14 uM", "inactive")
    .to_numpy()
    .astype(np.float32)
)
truth[:, 1:3] = np.clip(truth[:, 1:3], -150, 50)
y = truth[:, 0]
meta = compounds.select("compound_id").join(
    assignments, on="compound_id", how="left", validate="1:1"
)
folds = meta["fold"].fill_null(-1).to_numpy()
anchors = meta["is_anchor"].fill_null(False).to_numpy()
groups = meta["group"].cast(pl.String).fill_null("test").to_numpy()
f = featurizers.SimpleMoleculeMolGraphFeaturizer()
graphs = [f(Chem.MolFromSmiles(s)) for s in allc["standardized_smiles"]]


def graph(ix):
    b = BatchMolGraph([graphs[int(i)] for i in ix])
    b.to(device)
    return b


def evaluate(model, indices, means, stds):
    model.eval()
    outp = []
    with torch.inference_mode():
        for ix in np.array_split(indices, max(1, int(np.ceil(len(indices) / 128)))):
            outp.append(model(graph(ix)).cpu().numpy())
    return np.concatenate(outp) * stds + means


allrows = []
membership = []
stats = []
grads = []
for fold in range(5):
    eligible = np.isfinite(y[:6896]) & ((folds != fold) | anchors)
    outer = (folds == fold) & ~anchors
    rng = np.random.default_rng(42 + fold)
    fixed = set(groups[eligible & anchors])
    possible = sorted(set(groups[eligible]) - fixed)
    valgroups = set(rng.permutation(possible)[: max(1, int(len(possible) * 0.1))])
    inner_val = eligible & np.isin(groups, list(valgroups))
    inner_train = eligible & ~inner_val
    assert set(groups[inner_train]).isdisjoint(set(groups[inner_val]))
    ti = np.flatnonzero(inner_train)
    vi = np.flatnonzero(inner_val)
    if refit_epochs is not None:
        ti = np.flatnonzero(eligible)
        inner_train = eligible
        inner_val = np.zeros_like(eligible)
    ei = np.arange(6896, n)
    membership.append(
        meta.select("compound_id").with_columns(
            pl.lit(fold).alias("outer_fold"),
            pl.Series("training", inner_train),
            pl.Series("inner_validation", inner_val),
            pl.Series("outer_holdout", outer),
        )
    )
    lv = (
        labels.join(meta.select("compound_id"), on="compound_id", validate="1:1")
        .filter(pl.col("compound_id").is_in(compounds[vi]["compound_id"].implode()))
        .sort("compound_id")
    )
    vi = np.array([lookup[c] for c in lv["compound_id"]])
    for public_weight in cfg["public_weights"]:
        tag = f"p{int(public_weight * 100):02d}"
        directory = out / f"{tag}-fold{fold}"
        directory.mkdir()
        torch.manual_seed(42 + fold)
        rng = np.random.default_rng(42 + fold)
        trainall = np.concatenate([ti, ei])
        means = np.nanmean(truth[trainall], 0)
        stds = np.nanstd(truth[trainall], 0)
        stds = np.maximum(stds, 0.05)
        means[3] = 0
        stds[3] = 1
        z = (truth - means) / stds
        counts = np.isfinite(truth[trainall]).sum(0)
        mass = counts * np.array([1, public_weight, public_weight, public_weight])
        normalizer = float(mass.sum())
        (directory / "scaling.json").write_text(
            json.dumps(
                dict(
                    means=means.tolist(),
                    stds=stds.tolist(),
                    observed_label_counts=counts.tolist(),
                    weighted_label_mass=mass.tolist(),
                    public_mass_fraction=float(mass[1:].sum() / normalizer),
                    public_weight=public_weight,
                ),
                indent=2,
            )
        )
        model = MultitaskModel(state, 0.1, 4, arch, pretrained=True).to(device)
        optimizer = torch.optim.AdamW(
            [
                {"params": model.mp.parameters(), "lr": 1e-5},
                {"params": model.head.parameters(), "lr": 1e-4},
            ],
            weight_decay=1e-4,
        )
        best = float("inf")
        stale = 0
        best_epoch = 0
        epoch_budget = refit_epochs[(fold, tag)] if refit_epochs is not None else cfg["max_epochs"]
        for epoch in range(epoch_budget):
            model.train()
            steps = int(np.ceil(len(ti) / 32))
            db = np.array_split(rng.permutation(ti), steps)
            es = rng.choice(ei, size=min(cfg["public_samples_per_epoch"], len(ei)), replace=False)
            eb = np.array_split(es, steps)
            inclusion = len(ei) / len(es)
            total_d = total_e = 0.0
            for step, (didx, eidx) in enumerate(zip(db, eb, strict=True)):
                idx = np.concatenate([didx, eidx]) if public_weight else didx
                p = model(graph(idx))
                t = torch.tensor(np.nan_to_num(z[idx]), device=device)
                m = torch.tensor(np.isfinite(z[idx]), device=device)
                dloss = (((p[:, 0] - t[:, 0]) ** 2) * m[:, 0]).sum() * steps / normalizer
                error = (p[:, 1:3] - t[:, 1:3]) ** 2
                binary = torch.nn.functional.binary_cross_entropy_with_logits(
                    p[:, 3], t[:, 3], reduction="none"
                )
                eloss = (
                    ((error * m[:, 1:3]).sum() + (binary * m[:, 3]).sum())
                    * public_weight
                    * inclusion
                    * steps
                    / normalizer
                )
                if epoch == 0 and step == 0:
                    params = list(model.mp.parameters())
                    dg = torch.autograd.grad(dloss, params, retain_graph=True, allow_unused=True)
                    eg = torch.autograd.grad(eloss, params, retain_graph=True, allow_unused=True)

                    def norm(gg):
                        return float(
                            torch.sqrt(sum((g.detach() ** 2).sum() for g in gg if g is not None))
                        )

                    grads.append(
                        dict(
                            fold=fold,
                            method=tag,
                            direct_mp_gradient_norm=norm(dg),
                            public_mp_gradient_norm=norm(eg),
                            public_weight=public_weight,
                        )
                    )
                optimizer.zero_grad(set_to_none=True)
                loss = dloss + eloss
                assert torch.isfinite(loss)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 5, error_if_nonfinite=True)
                optimizer.step()
                total_d += float(dloss.detach()) / steps
                total_e += float(eloss.detach()) / steps
            pv = evaluate(model, vi, means, stds)[:, 0]
            score = (
                float(
                    scorer.rae_soft_threshold_absolute_error(
                        lv["value"].to_numpy(),
                        pv,
                        y_true_lower=lv["conf_low"].to_numpy(),
                        y_true_upper=lv["conf_high"].to_numpy(),
                    )
                )
                if refit_epochs is None
                else None
            )
            row = dict(
                method=tag,
                fold=fold,
                epoch=epoch + 1,
                training_direct_loss=total_d,
                training_public_loss=total_e,
                inner_ST_RAE=score,
                public_sample_n=len(es) if public_weight else 0,
                public_inverse_sampling_weight=inclusion if public_weight else 0,
            )
            with (directory / "history.jsonl").open("a") as stream:
                stream.write(json.dumps(row) + "\n")
            if refit_epochs is not None or score < best:
                best = score
                best_epoch = epoch + 1
                stale = 0
                best_state = copy.deepcopy(model.state_dict())
                torch.save(
                    dict(
                        state_dict=best_state,
                        means=means.tolist(),
                        stds=stds.tolist(),
                        config=cfg,
                        encoder_architecture=arch,
                        fold=fold,
                        public_weight=public_weight,
                        epoch=epoch + 1,
                    ),
                    directory / "best.pt",
                )
            else:
                stale += 1
            if (epoch + 1) % 5 == 0:
                print(json.dumps(row), flush=True)
            if refit_epochs is None and stale >= cfg["patience"]:
                break
        model.load_state_dict(best_state)
        model.eval()
        pred = evaluate(model, np.arange(6896), means, stds)
        with torch.inference_mode():
            latent = np.concatenate(
                [model.embed(graph(ix)).cpu().numpy() for ix in np.array_split(np.arange(6896), 54)]
            )
        assert np.isfinite(pred).all() and np.isfinite(latent).all()
        compounds.select("compound_id").with_columns(
            [pl.Series(f"emb_{j:03d}", latent[:, j]) for j in range(600)]
        ).write_parquet(directory / "embeddings.parquet")
        frame = labels.join(assignments, on="compound_id", validate="1:1").sort("compound_id")
        ids = np.array([lookup[c] for c in frame["compound_id"]])
        isval = (frame["fold"].to_numpy() == fold) & ~frame["is_anchor"].to_numpy()
        common = frame.filter(pl.Series(isval)).select(
            "compound_id",
            "endpoint",
            "fold",
            "is_neighbor",
            "conf_low",
            "conf_high",
            pl.col("value").alias("y_true"),
        )
        allrows.append(
            common.with_columns(
                pl.Series("prediction", pred[ids[isval], 0].astype(float)),
                pl.lit("d6_joint_" + tag + "_direct" + ("_refit" if refit_epochs else "")).alias(
                    "method"
                ),
                pl.lit(cfg["split_run"]).alias("split_id"),
            )
        )
        aux = (
            frame.select("compound_id")
            .join(
                nested.filter(pl.col("outer_fold") == fold),
                on="compound_id",
                maintain_order="left",
                validate="1:1",
            )
            .select([f"aux_{j}" for j in range(4)])
            .to_numpy()
        )
        xx = np.column_stack([latent[ids], aux])
        tree = lgb.LGBMRegressor(
            objective="regression_l1",
            n_estimators=300,
            num_leaves=15,
            min_child_samples=20,
            learning_rate=0.05,
            reg_lambda=1,
            n_jobs=2,
            random_state=42 + fold,
            verbosity=-1,
            deterministic=True,
            force_col_wise=True,
        ).fit(xx[~isval], frame["value"].to_numpy()[~isval])
        tree.booster_.save_model(str(directory / "readout.lgbm.txt"))
        allrows.append(
            common.with_columns(
                pl.Series("prediction", tree.predict(xx[isval])),
                pl.lit("d6_joint_" + tag + "_lgbm" + ("_refit" if refit_epochs else "")).alias(
                    "method"
                ),
                pl.lit(cfg["split_run"]).alias("split_id"),
            )
        )
        stats.append(
            dict(
                method=tag,
                fold=fold,
                epochs=epoch + 1,
                best_epoch=best_epoch,
                inner_ST_RAE=best,
                public_weight=public_weight,
                public_label_mass_fraction=float(mass[1:].sum() / normalizer),
            )
        )
        print(json.dumps(stats[-1]), flush=True)
        pl.concat(allrows).write_parquet(out / "joint-oof.parquet")
        pl.DataFrame(grads).write_csv(out / "encoder-gradient-audit.csv")
        (out / "progress.json").write_text(json.dumps(stats, indent=2))
pl.concat(membership).write_parquet(out / "memberships.parquet")
external.select("compound_id", "connectivity").write_parquet(out / "external-membership.parquet")
metrics = []
for (method,), g in pl.concat(allrows).group_by("method"):
    assert g.height == 1457 and g.select("compound_id", "endpoint", "fold").n_unique() == 1457
    metrics.append(
        dict(
            method=method,
            ST_RAE=float(
                scorer.rae_soft_threshold_absolute_error(
                    g["y_true"].to_numpy(),
                    g["prediction"].to_numpy(),
                    y_true_lower=g["conf_low"].to_numpy(),
                    y_true_upper=g["conf_high"].to_numpy(),
                )
            ),
        )
    )
pl.DataFrame(metrics).write_csv(out / "metrics.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            folds=stats,
            heads=[
                "D6_direct_pIC50",
                "AID1851_max_response",
                "AID1851_activity57.14uM",
                "AID1851_inactive",
            ],
            architecture=arch,
            outer_refit=refit_epochs is not None,
            interpretation=(
                "Joint encoder gradients; per-observed-label weighted standardized losses; "
                "sampled public inverse probability correction; "
                + (
                    "full outer refit at selected epochs"
                    if refit_epochs
                    else "selected inner checkpoint retained, no outer refit"
                )
            ),
        ),
        indent=2,
    )
)
