"""Matched multitask MLP versus constrained Hill observation-surrogate heads."""

import importlib.util
import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
import torch
from sklearn.model_selection import GroupShuffleSplit
from torch import nn

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


def read(i):
    return pl.read_parquet(path(i))


targets, assignments, emb, scores, labels, nested = [read(i) for i in range(6)]
spec = importlib.util.spec_from_file_location("official", path(6))
official = importlib.util.module_from_spec(spec)
spec.loader.exec_module(official)
torch.set_num_threads(2)
device = "cuda" if torch.cuda.is_available() else "cpu"
base = emb.join(scores, on="compound_id", validate="1:1").sort("compound_id")
ids = base["compound_id"].to_list()
lookup = {c: i for i, c in enumerate(ids)}
cols = [c for c in base.columns if c.startswith(("emb_", "score_"))]
x = base.select(cols).to_numpy().astype("float32")
assert x.shape == (6896, 647) and np.isfinite(x).all()
eps = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
kinds = ["direct", "tdi", "emax", "logfc"]
y = np.full((len(ids), 16), np.nan, dtype="float32")
lo = y.copy()
hi = y.copy()
concentration = np.full((len(ids), 4), np.nan, dtype="float32")
for r in targets.iter_rows(named=True):
    i = lookup[r["compound_id"]]
    j = 4 * kinds.index(r["kind"]) + eps.index(r["enzyme"])
    y[i, j] = r["value"]
    lo[i, j] = r["conf_low"] if r["conf_low"] is not None else r["value"]
    hi[i, j] = r["conf_high"] if r["conf_high"] is not None else r["value"]
    if r["kind"] == "logfc":
        concentration[i, j % 4] = r["concentration_M"]
for j in range(4):
    values = np.unique(concentration[np.isfinite(concentration[:, j]), j])
    assert len(values) == 1
    concentration[:, j] = values[0]
mask = np.isfinite(y)
yfill = np.nan_to_num(y)
lofill = np.nan_to_num(lo)
hifill = np.nan_to_num(hi)
a = base.select("compound_id").join(assignments, on="compound_id", how="left", validate="1:1")
folds = a["fold"].fill_null(-1).to_numpy()
anchors = a["is_anchor"].fill_null(False).to_numpy()
groups = a["group"].cast(pl.String).fill_null("unassigned").to_numpy()
available = (folds >= 0) & mask.any(axis=1)


class Model(nn.Module):
    def __init__(self, arm, mu, sd):
        super().__init__()
        self.arm = arm
        self.net = nn.Sequential(
            nn.Linear(647, 256),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Dropout(0.1),
        )
        self.head = nn.Linear(128, 16)
        self.register_buffer("mu", torch.tensor(mu))
        self.register_buffer("sd", torch.tensor(sd))
        self.logscale = nn.Parameter(torch.zeros(4))
        self.offset = nn.Parameter(torch.zeros(4))
        self.hraw = nn.Parameter(
            torch.full((4,), float(np.log(0.2 / 0.8))), requires_grad=arm == "hill_learned"
        )

    def forward(self, z, c):
        latent = self.net(z)
        p = self.head(latent) * self.sd + self.mu
        h = (
            0.25 + 3.75 * torch.sigmoid(self.hraw)
            if self.arm == "hill_learned"
            else torch.ones(4, device=z.device)
        )
        if self.arm != "independent":
            response = p[:, 8:12] * torch.sigmoid(np.log(10) * h * (p[:, 4:8] + torch.log10(c)))
            p = torch.cat(
                [p[:, :12], torch.nn.functional.softplus(self.logscale) * response + self.offset],
                dim=1,
            )
        return p, latent, h


def score(pred, indices):
    result = []
    for j in range(4):
        ii = indices[mask[indices, j]]
        result.append(
            float(
                official.rae_soft_threshold_absolute_error(
                    y[ii, j],
                    pred[mask[indices, j], j],
                    y_true_lower=lo[ii, j],
                    y_true_upper=hi[ii, j],
                )
            )
        )
    return float(np.mean(result))


logs = []
members = []
oofs = []
fits = []
for fold in range(5):
    outer = available & ((folds != fold) | anchors)
    val = available & (folds == fold) & ~anchors
    tridx = np.flatnonzero(outer)
    vaidx = np.flatnonzero(val)
    anchor_groups = set(groups[tridx[anchors[tridx]]])
    fixed_train = np.array([i for i in tridx if groups[i] in anchor_groups], dtype=int)
    nonanchor = np.array([i for i in tridx if groups[i] not in anchor_groups], dtype=int)
    it, iv = next(
        GroupShuffleSplit(n_splits=1, test_size=0.1, random_state=42 + fold).split(
            nonanchor, groups=groups[nonanchor]
        )
    )
    innertr = np.concatenate([nonanchor[it], fixed_train])
    innerva = nonanchor[iv]
    assert set(groups[innertr]).isdisjoint(set(groups[innerva]))
    members.append(
        pl.DataFrame(
            {
                "compound_id": ids,
                "outer_fold": [fold] * len(ids),
                "outer_train": outer,
                "outer_validation": val,
                "inner_train": np.isin(np.arange(len(ids)), innertr),
                "inner_validation": np.isin(np.arange(len(ids)), innerva),
            }
        )
    )
    assert not np.intersect1d(innertr, vaidx).size
    for arm in cfg["arms"]:
        selected = None
        final_model = None
        final_x = None
        for phase, trainidx, valididx in [("select", innertr, innerva), ("refit", tridx, vaidx)]:
            torch.manual_seed(42 + fold)
            np.random.seed(42 + fold)
            xm = x[trainidx].mean(0)
            xs = x[trainidx].std(0)
            xs = np.maximum(xs, 1e-6)
            mu = np.nanmean(y[trainidx], axis=0)
            sd = np.maximum(np.nanstd(y[trainidx], axis=0), 0.05)
            assert np.isfinite(mu).all() and np.isfinite(sd).all()
            denom = np.array(
                [
                    np.maximum(
                        np.maximum(lo[trainidx, j] - mu[j], 0)
                        + np.maximum(mu[j] - hi[trainidx, j], 0),
                        0,
                    )[mask[trainidx, j]].mean()
                    for j in range(4)
                ],
                dtype="float32",
            )
            denom = np.maximum(denom, 0.05)
            xx = torch.tensor((x - xm) / xs, device=device)
            yy = torch.tensor(yfill, device=device)
            ll = torch.tensor(lofill, device=device)
            hh = torch.tensor(hifill, device=device)
            mm = torch.tensor(mask, device=device)
            cc = torch.tensor(concentration, device=device)
            model = Model(arm, mu, sd).to(device)
            opt = torch.optim.AdamW(model.parameters(), lr=0.001, weight_decay=0.0001)
            dd = torch.tensor(denom, device=device)
            ss = torch.tensor(sd, device=device)
            best = float("inf")
            best_epoch = 1
            stale = 0
            for epoch in range(1, (cfg["max_epochs"] if phase == "select" else selected) + 1):
                model.train()
                losses = []
                for idx in np.array_split(
                    np.random.permutation(trainidx), int(np.ceil(len(trainidx) / 256))
                ):
                    p, _, h = model(xx[idx], cc[idx])
                    m = mm[idx]
                    direct = (
                        torch.relu(ll[idx, :4] - p[:, :4]) + torch.relu(p[:, :4] - hh[idx, :4])
                    ) / dd
                    loss = 2 * ((direct * m[:, :4]).sum(0) / m[:, :4].sum(0).clamp(min=1)).mean()
                    for start, weight in [(4, 0.5), (8, 0.25), (12, 0.5)]:
                        e = (
                            (p[:, start : start + 4] - yy[idx, start : start + 4])
                            / ss[start : start + 4]
                        ) ** 2
                        loss = (
                            loss
                            + weight
                            * (
                                (e * m[:, start : start + 4]).sum(0)
                                / m[:, start : start + 4].sum(0).clamp(min=1)
                            ).mean()
                        )
                    if arm == "hill_learned":
                        loss = loss + 0.01 * (torch.log(h) ** 2).mean()
                    opt.zero_grad()
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(model.parameters(), 5)
                    opt.step()
                    losses.append(float(loss.detach()))
                model.eval()
                if phase == "select":
                    with torch.no_grad():
                        vp = model(xx[valididx], cc[valididx])[0].cpu().numpy()
                    metric = score(vp, valididx)
                    if metric < best - 1e-6:
                        best = metric
                        best_epoch = epoch
                        stale = 0
                    else:
                        stale += 1
                else:
                    metric = None
                logs.append(
                    dict(
                        arm=arm,
                        fold=fold,
                        phase=phase,
                        epoch=epoch,
                        loss=float(np.mean(losses)),
                        inner_ST_RAE=metric,
                    )
                )
                if phase == "select" and stale >= cfg["patience"]:
                    break
            if phase == "select":
                selected = best_epoch
            else:
                final_model = model
                final_x = xx
                torch.save(
                    dict(
                        state_dict=model.state_dict(),
                        x_mean=xm,
                        x_std=xs,
                        head_mean=mu,
                        head_std=sd,
                        direct_denominator=denom,
                        arm=arm,
                        selected_epochs=selected,
                        feature_columns=cols,
                    ),
                    out / f"{arm}-fold{fold}.pt",
                )
                with torch.no_grad():
                    p, latent, h = model(xx, cc)
                    p = p.cpu().numpy()
                    latent = latent.cpu().numpy()
                assert np.isfinite(p).all() and np.isfinite(latent).all()
                pl.DataFrame(
                    {
                        "compound_id": ids,
                        "outer_fold": [fold] * len(ids),
                        **{f"emb_{j:03d}": latent[:, j] for j in range(128)},
                    }
                ).write_parquet(out / f"{arm}-fold{fold}-embeddings.parquet")
                fits.append(
                    dict(
                        arm=arm,
                        fold=fold,
                        epochs=selected,
                        inner_ST_RAE=best if phase == "select" else None,
                        hill=h.cpu().detach().tolist(),
                        outer_train=len(tridx),
                        outer_validation=len(vaidx),
                    )
                )
                for j, ep in enumerate(eps):
                    epname = ep + "_pIC50_direct_inhibition"
                    frame = (
                        labels.filter(pl.col("endpoint") == epname)
                        .join(assignments, on="compound_id", validate="1:1")
                        .filter((pl.col("fold") == fold) & ~pl.col("is_anchor"))
                        .sort("compound_id")
                    )
                    indices = np.array([lookup[c] for c in frame["compound_id"]])
                    assert set(indices) <= set(vaidx)
                    common = frame.select(
                        "compound_id",
                        "endpoint",
                        "fold",
                        "is_neighbor",
                        "conf_low",
                        "conf_high",
                        pl.col("value").alias("y_true"),
                    )
                    oofs.append(
                        common.with_columns(
                            pl.lit(arm + "_direct").alias("method"),
                            pl.lit(cfg["split_run"]).alias("split_id"),
                            pl.Series("prediction", p[indices, j].astype(float)),
                        )
                    )
                    head = (
                        labels.filter(pl.col("endpoint") == epname)
                        .join(assignments, on="compound_id", validate="1:1")
                        .sort("compound_id")
                    )
                    allidx = np.array([lookup[c] for c in head["compound_id"]])
                    isval = (head["fold"].to_numpy() == fold) & ~head["is_anchor"].to_numpy()
                    aux = (
                        head.select("compound_id")
                        .join(
                            nested.filter(pl.col("outer_fold") == fold),
                            on="compound_id",
                            maintain_order="left",
                            validate="1:1",
                        )
                        .select([f"aux_{k}" for k in range(4)])
                        .to_numpy()
                    )
                    rx = np.column_stack([latent[allidx], aux])
                    assert np.isfinite(rx).all()
                    lm = lgb.LGBMRegressor(
                        objective="regression_l1",
                        n_estimators=300,
                        num_leaves=15,
                        min_child_samples=20,
                        learning_rate=0.05,
                        reg_lambda=1,
                        n_jobs=3,
                        random_state=42 + fold,
                        verbosity=-1,
                        deterministic=True,
                        force_col_wise=True,
                    )
                    lm.fit(rx[~isval], head["value"].to_numpy()[~isval])
                    lm.booster_.save_model(str(out / f"{arm}-{ep}-fold{fold}.lgbm.txt"))
                    oofs.append(
                        common.with_columns(
                            pl.lit(arm + "_lgbm").alias("method"),
                            pl.lit(cfg["split_run"]).alias("split_id"),
                            pl.Series("prediction", lm.predict(rx[isval])),
                        )
                    )
                print(json.dumps(fits[-1]), flush=True)
        pl.DataFrame(logs).write_parquet(out / "training-history.parquet")
        pl.concat(oofs).write_parquet(out / "observation-oof.parquet")
pl.concat(members).write_parquet(out / "memberships.parquet")
res = pl.concat(oofs)
rows = []
for (method, ep), f in res.group_by("method", "endpoint"):
    assert (
        f.height
        == {"CYP1A2": 1372, "CYP2C9": 1255, "CYP2D6": 1457, "CYP3A4": 2300}[ep.split("_")[0]]
    )
    rows.append(
        dict(
            method=method,
            endpoint=ep,
            ST_RAE=float(
                official.rae_soft_threshold_absolute_error(
                    f["y_true"].to_numpy(),
                    f["prediction"].to_numpy(),
                    y_true_lower=f["conf_low"].to_numpy(),
                    y_true_upper=f["conf_high"].to_numpy(),
                )
            ),
        )
    )
pl.DataFrame(rows).write_csv(out / "metrics.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            fits=fits,
            features=cols,
            device=device,
            interpretation=(
                "Constrained surrogate observation model; affine log2fc map "
                "is learned and is not experimental calibration"
            ),
        ),
        indent=2,
    )
)
print(pl.DataFrame(rows))
