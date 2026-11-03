"""Continue MONROE pretraining on the challenge-disjoint external CYP corpus.

Every CYP assay with enough rows becomes a task with a regression output (exact
potency values, standardized per task) and a binary output (activity calls, including
inactive AC50 records). The corpus excludes challenge compounds, so the resulting
encoder can initialize every outer fold. A small compound holdout only monitors loss.
Weights are saved in the upstream `encoder.*` layout so `load_ckpt` can read them.
"""

import json
import math
import os
import pickle
import shutil
from pathlib import Path

import numpy as np
import polars as pl
import torch
from monroe.eval.embed import to_pyg
from monroe.model.ckpt import load_ckpt
from torch_geometric.data import Batch


def schedule(step, warmup, steps):
    """Linear warmup then cosine decay, as a learning-rate multiplier."""
    return min(1.0, (step + 1) / warmup) * 0.5 * (1 + math.cos(math.pi * min(1.0, step / steps)))


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    device = cfg["device"]
    torch.manual_seed(cfg["seed"])
    rng = np.random.default_rng(cfg["seed"])
    checkpoint = Path(os.environ["TMPDIR"]) / "checkpoint"
    checkpoint.mkdir(parents=True, exist_ok=True)
    for name, alias in cfg["checkpoint"].items():
        if not (checkpoint / name).exists():
            (checkpoint / name).symlink_to(inp / alias)
    with open(inp / cfg["graphs"], "rb") as handle:
        graphs = pickle.load(handle)

    obs = pl.read_parquet(inp / "observations.parquet").filter(
        pl.col("compound_id").is_in(list(graphs))
    )
    counts = obs.group_by("assay_id").len().filter(pl.col("len") >= cfg["min_rows"])
    tasks = sorted(counts["assay_id"].to_list())
    regression = (
        obs.filter(pl.col("assay_id").is_in(tasks) & pl.col("value").is_not_null())
        .filter(pl.col("relation") == "=")
        .group_by("compound_id", "assay_id")
        .agg(pl.col("value").mean())
    )
    binary = (
        obs.filter(pl.col("assay_id").is_in(tasks) & pl.col("binary_label").is_not_null())
        .group_by("compound_id", "assay_id")
        .agg(pl.col("binary_label").mean().round(0))
    )
    ids = sorted(set(regression["compound_id"]) | set(binary["compound_id"]))
    position = {cid: i for i, cid in enumerate(ids)}
    task_index = {task: j for j, task in enumerate(tasks)}
    y_reg = np.full((len(ids), len(tasks)), np.nan, dtype=np.float32)
    y_bin = np.full((len(ids), len(tasks)), np.nan, dtype=np.float32)
    for frame, matrix, column in ((regression, y_reg, "value"), (binary, y_bin, "binary_label")):
        rows = [position[c] for c in frame["compound_id"]]
        cols = [task_index[t] for t in frame["assay_id"]]
        matrix[rows, cols] = frame[column].to_numpy()
    mean, std = np.nanmean(y_reg, axis=0), np.nanstd(y_reg, axis=0)
    mean, std = np.nan_to_num(mean), np.where(np.isfinite(std) & (std > 0), std, 1.0)
    y_reg = (y_reg - mean) / std
    order = rng.permutation(len(ids))
    holdout = set(order[: int(cfg["holdout_fraction"] * len(ids))].tolist())
    train_rows = [i for i in range(len(ids)) if i not in holdout]
    check_rows = sorted(holdout)

    encoder = load_ckpt(str(checkpoint)).to(device)
    head = torch.nn.Sequential(
        torch.nn.Dropout(cfg["dropout"]),
        torch.nn.Linear(cfg["embedding_dim"], cfg["head_hidden"]),
        torch.nn.GELU(),
        torch.nn.Linear(cfg["head_hidden"], 2 * len(tasks)),
    ).to(device)
    optimizer = torch.optim.AdamW(
        [
            {"params": encoder.parameters(), "lr": cfg["lr_encoder"]},
            {"params": head.parameters(), "lr": cfg["lr_head"]},
        ],
        weight_decay=cfg["weight_decay"],
    )
    steps = cfg["epochs"] * math.ceil(len(train_rows) / cfg["batch_size"])
    warmup = max(1, int(cfg["warmup_fraction"] * steps))
    scheduler = torch.optim.lr_scheduler.LambdaLR(
        optimizer, lambda s, w=warmup, n=steps: schedule(s, w, n)
    )
    t_reg, t_bin = torch.tensor(y_reg, device=device), torch.tensor(y_bin, device=device)

    def loss_for(rows):
        batch = Batch.from_data_list([to_pyg(graphs[ids[r]]) for r in rows]).to(device)
        index = torch.tensor(rows, device=device)
        with torch.autocast("cuda", dtype=torch.bfloat16):
            output = head(encoder(batch)[0]).float()
        reg, logit = output[:, : len(tasks)], output[:, len(tasks) :]
        yr, yb = t_reg[index], t_bin[index]
        mr, mb = torch.isfinite(yr), torch.isfinite(yb)
        squared = torch.where(mr, (reg - torch.nan_to_num(yr)) ** 2, 0.0)
        bce = torch.nn.functional.binary_cross_entropy_with_logits(
            logit, torch.nan_to_num(yb), reduction="none"
        )
        bce = torch.where(mb, bce, 0.0)
        return (squared.sum() + cfg["binary_weight"] * bce.sum()) / (mr.sum() + mb.sum()).clamp(
            min=1
        )

    log = []
    for epoch in range(cfg["epochs"]):
        encoder.train()
        head.train()
        shuffled = rng.permutation(train_rows)
        losses = []
        for start in range(0, len(shuffled), cfg["batch_size"]):
            loss = loss_for(shuffled[start : start + cfg["batch_size"]].tolist())
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(
                list(encoder.parameters()) + list(head.parameters()), cfg["grad_clip"]
            )
            optimizer.step()
            scheduler.step()
            losses.append(loss.item())
        encoder.eval()
        head.eval()
        with torch.no_grad():
            check = [
                loss_for(check_rows[s : s + cfg["batch_size"]]).item()
                for s in range(0, len(check_rows), cfg["batch_size"])
            ]
        log.append(
            {"epoch": epoch, "train": float(np.mean(losses)), "holdout": float(np.mean(check))}
        )
        print(json.dumps(log[-1]), flush=True)

    state = {f"encoder.{k}": v.detach().cpu() for k, v in encoder.state_dict().items()}
    torch.save(state, out / "weights.pt")
    shutil.copyfile(inp / cfg["checkpoint"]["config.json"], out / "config.json")
    reloaded = load_ckpt(str(out))
    for key, value in reloaded.state_dict().items():
        assert torch.equal(value, state[f"encoder.{key}"].to(value.dtype))
    summary = {
        "tasks": tasks,
        "compounds": len(ids),
        "holdout": len(check_rows),
        "regression_labels": int(np.isfinite(y_reg).sum()),
        "binary_labels": int(np.isfinite(y_bin).sum()),
        "log": log,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
