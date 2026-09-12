"""Execute the published Workbench Chemprop trainer on immutable author data.

The original source functions construct and train the native Chemprop MPNN.
This adapter supplies Polars-backed frames, canonical evaluation identities,
and durable training instrumentation; it does not replace the loss or optimizer.
"""

import ast
import gc
import importlib.util
import json
import os
import sys
import time
import types
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import models
from lightning import pytorch as lightning


class Frame:
    """Minimal read-only frame surface used by original trainer, without pandas."""

    def __init__(self, frame):
        self.frame = frame

    def __getitem__(self, column):
        return self.frame[column].to_numpy()

    def __len__(self):
        return self.frame.height


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class Trace(lightning.Callback):
    def __init__(self, output, fold):
        self.output, self.fold = output, fold
        self.start = time.monotonic()
        self.epoch_start = self.start

    def on_train_epoch_start(self, trainer, model):
        self.epoch_start = time.monotonic()

    def on_train_epoch_end(self, trainer, model):
        metrics = {
            key: float(value.detach().cpu())
            for key, value in trainer.callback_metrics.items()
            if isinstance(value, torch.Tensor) and value.numel() == 1
        }
        record = dict(
            fold=self.fold,
            epoch=trainer.current_epoch,
            seconds=time.monotonic() - self.epoch_start,
            elapsed_seconds=time.monotonic() - self.start,
            global_step=trainer.global_step,
            training_batches=trainer.num_training_batches,
            lr=[group["lr"] for group in trainer.optimizers[0].param_groups],
            gpu_peak_allocated=torch.cuda.max_memory_allocated()
            if torch.cuda.is_available()
            else 0,
            **metrics,
        )
        with (self.output / "history.jsonl").open("a") as stream:
            stream.write(json.dumps(record) + "\n")
        if trainer.current_epoch < 3 or trainer.current_epoch % 10 == 0:
            print(json.dumps(record), flush=True)


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    for package in ["workbench", "workbench.endpoints", "workbench.training"]:
        module = types.ModuleType(package)
        module.__path__ = []
        sys.modules[package] = module
    load_module("workbench.endpoints.chemprop_utils", inp / cfg["source_helpers"])
    core = load_module("workbench.training.chemprop_core", inp / cfg["source_core"])
    splits = load_module("workbench.training.splits", inp / cfg["source_splits"])
    syntax = ast.parse((inp / cfg["source_template"]).read_text())
    hp = next(
        ast.literal_eval(node.value)
        for node in syntax.body
        if isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "DEFAULT_HYPERPARAMETERS" for t in node.targets)
    )
    expected = dict(
        hidden_dim=700,
        depth=5,
        ffn_hidden_dim=300,
        ffn_num_layers=2,
        dropout=0.1,
        criterion="mae",
        max_epochs=400,
        patience=40,
        batch_size=64,
        warmup_epochs=2,
        init_lr=1e-4,
        max_lr=1e-3,
        final_lr=1e-4,
        seed=42,
        from_foundation=None,
        n_folds=5,
    )
    assert all(hp[key] == value for key, value in expected.items())
    torch.set_num_threads(cfg.get("num_threads", 2))
    torch.set_float32_matmul_precision("medium")
    table = pl.read_parquet(inp / cfg["table_file"])
    table = table.with_columns(pl.col("compound_id").alias("row_id"))
    heads = cfg["heads"]
    # Source filters all-head-missing rows before scaffold grouping.
    table = table.filter(pl.any_horizontal([pl.col(head).is_finite() for head in heads]))
    assert table["row_id"].n_unique() == table.height
    assert table["smiles"].null_count() == 0
    y = table.select(heads).to_numpy()
    smiles = table["smiles"].to_list()
    validated, indices = core.create_molecule_datapoints(smiles, y)
    assert len(indices) == table.height, (
        "Author table contains invalid SMILES; explicit curation needed"
    )
    del validated
    n_primary = cfg["n_primary"]
    counts = np.isfinite(y).sum(axis=0)
    assert (counts > 0).all()
    if cfg["recipe"] in ["union18", "aux8"]:
        inv = 1.0 / counts[:n_primary].astype(np.float64)
        primary_weights = inv / inv.mean()
        task_weights = np.r_[
            primary_weights, np.full(len(heads) - n_primary, 0.3 * primary_weights.mean())
        ]
    elif cfg["recipe"] == "isoform6":
        task_weights = np.array([1.0] + [0.3] * (len(heads) - 1))
    else:
        assert cfg["recipe"] == "single1" and len(heads) == 1
        task_weights = np.array([1.0])
    (out / "effective-config.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                hyperparameters=hp,
                heads=heads,
                global_head_counts=counts.tolist(),
                task_weights=task_weights.tolist(),
                target_weight_basis="Whole input table counts before CV, as published",
            ),
            indent=2,
        )
    )
    groups = splits.get_scaffold_groups(smiles)
    native = np.full(table.height, -1, dtype=int)
    for fold, (_, val) in enumerate(splits.group_folds(groups, 5, 42)):
        native[val] = fold
    assigned = pl.read_parquet(inp / cfg["assignments_file"])
    connection_map = pl.read_parquet(inp / cfg["connectivity_folds_file"]).filter(
        pl.col("fold").is_not_null()
    )
    conflicts = connection_map.group_by("connectivity").agg(pl.col("fold").n_unique())
    assert conflicts["fold"].max() == 1, "A connectivity has conflicting parent folds"
    connection_map = connection_map.group_by("connectivity").agg(
        pl.col("fold").first(), pl.col("is_anchor").any()
    )
    canonical = table.select("row_id", "compound_id", "connectivity").join(
        connection_map, on="connectivity", how="left", validate="m:1", maintain_order="left"
    )
    parent_fold = canonical["fold"].fill_null(-1).to_numpy()
    anchors = canonical["is_anchor"].fill_null(False).to_numpy()
    fold_ids = native.copy()
    if cfg.get("cv_policy", "parent_adapted") == "parent_adapted":
        fold_ids[parent_fold >= 0] = parent_fold[parent_fold >= 0]
        fold_ids[anchors] = -1
    else:
        assert cfg["cv_policy"] == "source_scaffold"
    # The prepared data maps every challenge connectivity to its canonical identity.
    # A source row may carry challenge and public labels; this single row is excluded as a unit.
    table.select("row_id", "compound_id", "smiles").with_columns(
        pl.Series("native_scaffold_group", groups),
        pl.Series("native_fold", native),
        pl.Series("parent_fold", parent_fold),
        pl.Series("is_anchor", anchors),
        pl.Series("effective_fold", fold_ids),
    ).write_parquet(out / "assignments.parquet")
    labels = pl.read_parquet(inp / cfg["labels_file"])
    assert set(cfg["scored_endpoints"]).issubset(set(labels["endpoint"]))
    test = pl.read_parquet(inp / cfg["test_file"])
    test_mapping = pl.read_parquet(inp / cfg["test_mapping_file"])
    test = test.select("compound_id", "smiles", "connectivity").join(
        test_mapping.select("compound_id", pl.col("molecule_name").alias("Molecule_Name")),
        on="compound_id",
        how="inner",
        validate="1:1",
        maintain_order="left",
    )
    assert test.height == 750 and test["compound_id"].n_unique() == 750
    overlap = table.filter(pl.col("connectivity").is_in(test["connectivity"].implode()))
    overlap.write_parquet(out / "public-blind-connectivity-overlap.parquet")
    oofs, tests, summary = [], [], []
    for fold in range(5):
        directory = out / f"fold{fold}"
        directory.mkdir()
        train_idx, val_idx = np.flatnonzero(fold_ids != fold), np.flatnonzero(fold_ids == fold)
        assert len(val_idx) and not set(train_idx) & set(val_idx)
        table.select("row_id", "compound_id").with_columns(
            pl.lit(fold).alias("outer_fold"),
            pl.Series("training", fold_ids != fold),
            pl.Series("validation", fold_ids == fold),
        ).write_parquet(directory / "membership.parquet")
        train_counts, val_counts = np.isfinite(y[train_idx]).sum(0), np.isfinite(y[val_idx]).sum(0)
        assert (train_counts > 0).all()
        trace = Trace(directory, fold)

        def trainer_factory(*args, _trace=trace, **kwargs):
            kwargs["callbacks"] = list(kwargs.get("callbacks", [])) + [_trace]
            return lightning.Trainer(*args, **kwargs)

        core.pl = types.SimpleNamespace(
            seed_everything=lightning.seed_everything,
            callbacks=lightning.callbacks,
            Trainer=trainer_factory,
        )
        spec = core.FoldSpec(
            hyperparameters=hp,
            task="regression",
            model_type="uq_regressor",
            smiles_column="smiles",
            n_targets=len(heads),
            task_weights=task_weights.astype(np.float32) if len(heads) > 1 else None,
            num_workers=cfg.get("num_workers", 0),
            concurrent=True,
            enable_progress_bar=False,
            verbose=True,
        )
        model, val_predictions = core.train_chemprop_fold(
            spec,
            Frame(table[train_idx]),
            Frame(table[val_idx]),
            y[train_idx],
            y[val_idx],
            fold_idx=fold,
            checkpoint_dir=str(directory),
        )
        assert (
            val_predictions.shape == (len(val_idx), len(heads))
            and np.isfinite(val_predictions).all()
        )
        models.save_model(directory / "model.pt", model, output_columns=heads)
        checkpoint_paths = list(directory.glob("best_*.ckpt"))
        assert len(checkpoint_paths) == 1
        checkpoint = torch.load(checkpoint_paths[0], map_location="cpu", weights_only=False)
        best_epoch = int(checkpoint["epoch"])
        pl.DataFrame(val_predictions, schema=heads).insert_column(
            0, table[val_idx]["row_id"]
        ).write_parquet(directory / "all-head-validation.parquet")
        test_predictions = core.predict_chemprop_frame(model, spec, Frame(test))
        assert test_predictions.shape == (750, len(heads)) and np.isfinite(test_predictions).all()
        pl.DataFrame(test_predictions, schema=heads).insert_column(
            0, test["compound_id"]
        ).write_parquet(directory / "all-head-test.parquet")
        for j, endpoint in enumerate(cfg["scored_endpoints"]):
            source = (
                table[val_idx]
                .select("compound_id")
                .with_columns(pl.Series("prediction", val_predictions[:, j].astype(float)))
                .filter(pl.col("compound_id").is_not_null())
            )
            evaluated = labels.filter(pl.col("endpoint") == endpoint).join(
                source, on="compound_id", how="inner", validate="1:1"
            )
            if cfg.get("cv_policy", "parent_adapted") == "parent_adapted":
                expected = (
                    labels.filter(pl.col("endpoint") == endpoint)
                    .join(assigned, on="compound_id")
                    .filter((pl.col("fold") == fold) & ~pl.col("is_anchor"))
                )
                assert set(evaluated["compound_id"]) == set(expected["compound_id"])
            oofs.append(
                evaluated.select(
                    "compound_id",
                    "endpoint",
                    "conf_low",
                    "conf_high",
                    pl.col("value").alias("y_true"),
                    "prediction",
                )
                .join(
                    assigned.select("compound_id", "is_neighbor"),
                    on="compound_id",
                    how="left",
                    validate="m:1",
                )
                .with_columns(
                    pl.lit(fold).alias("fold"),
                    pl.lit(cfg["split_run"]).alias("split_id"),
                    pl.lit(cfg["recipe"]).alias("method"),
                )
            )
            tests.append(
                test.select("compound_id", "Molecule_Name", "smiles").with_columns(
                    pl.lit(endpoint).alias("endpoint"),
                    pl.lit(fold).alias("fold"),
                    pl.Series("prediction", test_predictions[:, j].astype(float)),
                    pl.lit(cfg["recipe"]).alias("method"),
                )
            )
        summary.append(
            dict(
                fold=fold,
                train_rows=len(train_idx),
                validation_rows=len(val_idx),
                train_head_counts=train_counts.tolist(),
                val_head_counts=val_counts.tolist(),
                best_epoch=best_epoch,
                checkpoint_file=checkpoint_paths[0].name,
                seconds=time.monotonic() - trace.start,
            )
        )
        (out / "progress.json").write_text(json.dumps(summary, indent=2))
        pl.concat(oofs).write_parquet(out / "scp-oof.parquet")
        pl.concat(tests).write_parquet(out / "fold-test-predictions.parquet")
        del model, checkpoint
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        print(json.dumps(summary[-1]), flush=True)
    mean = (
        pl.concat(tests)
        .group_by("compound_id", "Molecule_Name", "smiles", "endpoint", "method")
        .agg(pl.col("prediction").mean(), pl.len().alias("n_folds"))
    )
    assert mean["n_folds"].min() == mean["n_folds"].max() == 5
    expected_oof = (
        (6384 if n_primary == 4 else 1457)
        if cfg.get("cv_policy", "parent_adapted") == "parent_adapted"
        else (6525 if n_primary == 4 else 1493)
    )
    assert pl.concat(oofs).height == expected_oof
    mean.write_parquet(out / "test-predictions.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                recipe=cfg["recipe"],
                source_commit=cfg["source_commit"],
                folds=summary,
                oof_rows=pl.concat(oofs).height,
                test_rows=mean.height,
                training="Original published train_chemprop_fold, native ChempropMAE and defaults",
                cv_policy=cfg.get("cv_policy", "parent_adapted"),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
