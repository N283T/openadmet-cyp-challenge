"""Fold-wise multi-task fine-tuning of the MONROE encoder on challenge regression tasks.

Main tasks are the four direct-inhibition pIC50 endpoints; optional auxiliary tasks
are other challenge measurements (TDI-condition pIC50, Emax) and the public primary
single-concentration log2fc. Missing labels are masked. Each outer parent fold trains
on compounds outside that fold for a fixed number of epochs (no validation-based
stopping), sampling one cached conformer per molecule per step. Validation and blind
compounds are predicted as the mean over all cached conformers; the blind member is
the five-fold mean. Fold-specific fine-tuned embeddings are also read out with TabPFN.
"""

import importlib.util
import json
import math
import os
import pickle
from pathlib import Path

import numpy as np
import polars as pl
import torch
from monroe.eval.embed import to_pyg
from monroe.model.ckpt import load_ckpt
from tabpfn import TabPFNRegressor
from tabpfn.constants import ModelVersion
from torch_geometric.data import Batch


def schedule(step, warmup, steps):
    """Linear warmup then cosine decay, as a learning-rate multiplier."""
    return min(1.0, (step + 1) / warmup) * 0.5 * (1 + math.cos(math.pi * min(1.0, step / steps)))


class Head(torch.nn.Module):
    def __init__(self, dim, tasks, hidden, dropout):
        super().__init__()
        self.net = torch.nn.Sequential(
            torch.nn.Dropout(dropout),
            torch.nn.Linear(dim, hidden),
            torch.nn.GELU(),
            torch.nn.Dropout(dropout),
            torch.nn.Linear(hidden, tasks),
        )

    def forward(self, x):
        return self.net(x)


def label_matrix(inp, cfg, ids):
    """Wide dev-compound task matrix keyed by compound ID; NaN marks missing labels."""
    measurements = pl.read_parquet(inp / "measurements.parquet")
    frame = pl.DataFrame({"compound_id": ids})
    for task in cfg["tasks"]:
        if task["source"] == "measurement":
            values = measurements.filter(pl.col("endpoint") == task["endpoint"]).select(
                "compound_id", pl.col("value").alias(task["name"])
            )
        else:
            aliases = (
                pl.read_parquet(inp / "aliases.parquet")
                .filter(pl.col("role") != "test")
                .select(pl.col("molecule_name").alias("Molecule_Name"), "compound_id")
                .unique()
            )
            values = (
                pl.read_csv(inp / "primary.csv")
                .filter(pl.col("enzyme") == task["enzyme"])
                .join(aliases, on="Molecule_Name", how="inner")
                .group_by("compound_id")
                .agg(pl.col("log2fc_estimate").mean().alias(task["name"]))
            )
        assert values["compound_id"].n_unique() == values.height
        frame = frame.join(values, on="compound_id", how="left", maintain_order="left")
    return frame


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    device = cfg["device"]
    torch.manual_seed(cfg["seed"])
    rng = np.random.default_rng(cfg["seed"])
    spec = importlib.util.spec_from_file_location("official", inp / "custom_scoring_functions.py")
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)

    def st_rae(frame, prediction):
        return float(
            scorer.rae_soft_threshold_absolute_error(
                frame["y_true"].to_numpy(),
                prediction,
                y_true_upper=frame["conf_high"].to_numpy(),
                y_true_lower=frame["conf_low"].to_numpy(),
            )
        )

    checkpoint = Path(os.environ["TMPDIR"]) / "checkpoint"
    checkpoint.mkdir(parents=True, exist_ok=True)
    for name, alias in cfg["checkpoint"].items():
        if not (checkpoint / name).exists():
            (checkpoint / name).symlink_to(inp / alias)
    graphs = []
    for name in cfg["graph_files"]:
        with open(inp / name, "rb") as handle:
            graphs.append(pickle.load(handle))
    assignments = pl.read_parquet(inp / "parent-assignments.parquet").sort("compound_id")
    blind_ids = pl.read_parquet(inp / "test-mapping.parquet")["compound_id"].sort().to_list()
    dev_ids = assignments["compound_id"].to_list()
    all_ids = dev_ids + blind_ids
    assert all(cid in graphs[0] for cid in all_ids)
    labels = label_matrix(inp, cfg, dev_ids)
    task_names = [task["name"] for task in cfg["tasks"]]
    task_weights = torch.tensor([task.get("weight", 1.0) for task in cfg["tasks"]], device=device)
    y_all = labels.select(task_names).to_numpy().astype(np.float32)
    dev_fold = assignments["fold"].to_numpy()
    (out / "task-counts.json").write_text(
        json.dumps({n: int(np.isfinite(y_all[:, i]).sum()) for i, n in enumerate(task_names)})
    )

    controls = {}
    for endpoint, control in cfg["controls"].items():
        frame = pl.read_parquet(inp / control["file"])
        for column, value in control["filter"].items():
            frame = frame.filter(pl.col(column) == value)
        frame = frame.filter(pl.col("endpoint") == endpoint).sort("compound_id")
        assert frame.height == control["rows"]
        controls[endpoint] = frame
    main_index = {task["endpoint"]: i for i, task in enumerate(cfg["tasks"]) if task.get("main")}
    assert set(main_index) == set(controls)

    def batches(ids, size, seed_choice):
        for start in range(0, len(ids), size):
            chunk = ids[start : start + size]
            data = [to_pyg(graphs[seed_choice(i)][cid]) for i, cid in enumerate(chunk)]
            yield chunk, Batch.from_data_list(data).to(device)

    def encode_all(encoder, ids):
        """Mean graph embedding over every cached conformer."""
        encoder.eval()
        total = None
        with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16):
            for conformer in range(len(graphs)):
                parts = [
                    encoder(batch)[0].float().cpu()
                    for _, batch in batches(ids, cfg["eval_batch_size"], lambda _, c=conformer: c)
                ]
                stacked = torch.cat(parts)
                total = stacked if total is None else total + stacked
        return (total / len(graphs)).numpy()

    oof_rows, blind_rows, records, readout_rows = [], [], [], []
    blind_head = {endpoint: [] for endpoint in controls}
    for fold in range(5):
        encoder = load_ckpt(str(checkpoint)).to(device)
        train_mask = (dev_fold != fold) & np.isfinite(y_all).any(axis=1)
        train_ids = [dev_ids[i] for i in np.flatnonzero(train_mask)]
        y_train = y_all[train_mask]
        mean, std = np.nanmean(y_train, axis=0), np.nanstd(y_train, axis=0)
        head = Head(cfg["embedding_dim"], len(task_names), cfg["head_hidden"], cfg["dropout"]).to(
            device
        )
        optimizer = torch.optim.AdamW(
            [
                {"params": encoder.parameters(), "lr": cfg["lr_encoder"]},
                {"params": head.parameters(), "lr": cfg["lr_head"]},
            ],
            weight_decay=cfg["weight_decay"],
        )
        steps = cfg["epochs"] * math.ceil(len(train_ids) / cfg["batch_size"])
        warmup = max(1, int(cfg["warmup_fraction"] * steps))
        scheduler = torch.optim.lr_scheduler.LambdaLR(
            optimizer, lambda s, w=warmup, n=steps: schedule(s, w, n)
        )
        target = torch.tensor((y_train - mean) / std, device=device)
        position = {cid: i for i, cid in enumerate(train_ids)}
        for epoch in range(cfg["epochs"]):
            encoder.train()
            head.train()
            order = rng.permutation(len(train_ids))
            shuffled = [train_ids[i] for i in order]
            draws = rng.integers(0, len(graphs), len(shuffled))
            losses = []
            for start in range(0, len(shuffled), cfg["batch_size"]):
                chunk = shuffled[start : start + cfg["batch_size"]]
                data = [to_pyg(graphs[draws[start + i]][cid]) for i, cid in enumerate(chunk)]
                batch = Batch.from_data_list(data).to(device)
                rows = torch.tensor([position[cid] for cid in chunk], device=device)
                y = target[rows]
                mask = torch.isfinite(y)
                with torch.autocast("cuda", dtype=torch.bfloat16):
                    prediction = head(encoder(batch)[0]).float()
                error = torch.where(mask, (prediction - torch.nan_to_num(y)) ** 2, 0.0)
                per_task = error.sum(0) / mask.sum(0).clamp(min=1)
                present = (mask.sum(0) > 0).float()
                loss = (per_task * task_weights * present).sum() / (task_weights * present).sum()
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(
                    list(encoder.parameters()) + list(head.parameters()), cfg["grad_clip"]
                )
                optimizer.step()
                scheduler.step()
                losses.append(loss.item())
            print(
                json.dumps({"fold": fold, "epoch": epoch, "loss": float(np.mean(losses))}),
                flush=True,
            )

        embeddings = encode_all(encoder, all_ids)
        np.save(out / f"embeddings-fold{fold}.npy", embeddings.astype(np.float32))
        head.eval()
        with torch.no_grad():
            standardized = head(torch.tensor(embeddings, device=device)).cpu().numpy()
        predictions = standardized * std + mean
        index = {cid: i for i, cid in enumerate(all_ids)}
        for endpoint, task in main_index.items():
            control = controls[endpoint]
            valid = control.filter(pl.col("fold") == fold)
            rows = [index[cid] for cid in valid["compound_id"]]
            oof_rows.append(
                valid.select(
                    "compound_id", "endpoint", "fold", "y_true", "conf_low", "conf_high"
                ).with_columns(
                    pl.lit("monroe_ft_head").alias("method"),
                    pl.Series("prediction", predictions[rows, task]),
                )
            )
            blind_head[endpoint].append(predictions[[index[c] for c in blind_ids], task])

            # TabPFN readout on this fold's fine-tuned embedding.
            labelled = (
                pl.read_parquet(inp / "measurements.parquet")
                .filter(pl.col("endpoint") == endpoint)
                .join(assignments.select("compound_id", "fold"), on="compound_id")
                .filter(pl.col("fold") != fold)
            )
            x_train = embeddings[[index[c] for c in labelled["compound_id"]]]
            model = TabPFNRegressor.create_default_for_version(
                ModelVersion.V2_6,
                n_estimators=cfg["tabpfn_estimators"],
                auto_scale_n_estimators=False,
                model_path=str(inp / "tabpfn.ckpt"),
                device=device,
                random_state=cfg["seed"] + fold,
                n_preprocessing_jobs=1,
                memory_saving_mode="auto",
                show_progress_bar=False,
                ignore_pretraining_limits=True,
                fit_mode="fit_preprocessors",
            )
            model.fit(x_train, labelled["value"].to_numpy())
            readout = np.asarray(model.predict(embeddings[rows]), dtype=float)
            readout_blind = np.asarray(
                model.predict(embeddings[[index[c] for c in blind_ids]]), dtype=float
            )
            readout_rows.append(
                valid.select(
                    "compound_id", "endpoint", "fold", "y_true", "conf_low", "conf_high"
                ).with_columns(
                    pl.lit("monroe_ft_tabpfn").alias("method"),
                    pl.Series("prediction", readout),
                )
            )
            blind_rows.append(
                pl.DataFrame({"compound_id": blind_ids}).with_columns(
                    pl.lit(endpoint).alias("endpoint"),
                    pl.lit("monroe_ft_tabpfn").alias("method"),
                    pl.lit(fold).alias("fold"),
                    pl.Series("prediction", readout_blind),
                )
            )
            del model
        # Persist this fold before the next one so a crash cannot discard it.
        pl.concat([f for f in oof_rows + readout_rows if f["fold"][0] == fold]).write_parquet(
            out / f"oof-fold{fold}.parquet"
        )
        pl.concat(
            [f for f in blind_rows if f["fold"][0] == fold]
            + [
                pl.DataFrame({"compound_id": blind_ids}).with_columns(
                    pl.lit(endpoint).alias("endpoint"),
                    pl.lit("monroe_ft_head").alias("method"),
                    pl.lit(fold).alias("fold"),
                    pl.Series("prediction", values[fold]),
                )
                for endpoint, values in blind_head.items()
            ]
        ).write_parquet(out / f"blind-fold{fold}.parquet")
        del encoder, head, optimizer
        torch.cuda.empty_cache()

    for endpoint, values in blind_head.items():
        for fold, prediction in enumerate(values):
            blind_rows.append(
                pl.DataFrame({"compound_id": blind_ids}).with_columns(
                    pl.lit(endpoint).alias("endpoint"),
                    pl.lit("monroe_ft_head").alias("method"),
                    pl.lit(fold).alias("fold"),
                    pl.Series("prediction", prediction),
                )
            )
    oof = pl.concat(oof_rows + readout_rows)
    oof.write_parquet(out / "oof.parquet")
    pl.concat(blind_rows).write_parquet(out / "blind-fold-predictions.parquet")

    for endpoint, control in controls.items():
        incumbent = control["prediction"].to_numpy()
        folds = control["fold"].to_numpy()
        for method in ("monroe_ft_head", "monroe_ft_tabpfn"):
            member = (
                control.select("compound_id")
                .join(
                    oof.filter((pl.col("endpoint") == endpoint) & (pl.col("method") == method)),
                    on="compound_id",
                    how="left",
                )["prediction"]
                .to_numpy()
            )
            record = {
                "endpoint": endpoint,
                "method": method,
                "ST_RAE": st_rae(control, member),
                "incumbent": st_rae(control, incumbent),
                "spearman": float(
                    pl.DataFrame({"a": member, "b": control["y_true"]})
                    .select(pl.corr("a", "b", method="spearman"))
                    .item()
                ),
                "residual_corr": float(
                    np.corrcoef(
                        control["y_true"].to_numpy() - member,
                        control["y_true"].to_numpy() - incumbent,
                    )[0, 1]
                ),
            }
            for weight in cfg["blend_weights"]:
                blended = (1 - weight) * incumbent + weight * member
                record[f"blend_{weight}"] = st_rae(control, blended)
                record[f"blend_{weight}_fold_wins"] = int(
                    sum(
                        st_rae(control.filter(pl.Series(folds == f)), blended[folds == f])
                        < st_rae(control.filter(pl.Series(folds == f)), incumbent[folds == f])
                        for f in range(5)
                    )
                )
            records.append(record)
            print(json.dumps(record), flush=True)
    pl.DataFrame(records).write_parquet(out / "metrics.parquet")
    (out / "summary.json").write_text(
        json.dumps({"config": cfg, "metrics": records}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
