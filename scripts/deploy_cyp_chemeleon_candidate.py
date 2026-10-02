"""Fit the validated equal-seed five-bag policy, persist weights, and score blind molecules."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from evaluate_cyp_d6_zero_start import fit_model, folds, predict
from validate_cyp_d6_bag_policy import select_cut


def deploy(table, fp, x, y, tdi, shift, weights, subtype, tasks, development, blind, cfg):
    torch.set_num_threads(cfg["cpu_threads"])
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    assert not np.isfinite(tasks[blind]).any()
    strategy, kind = cfg["deployment_split"], cfg["deployment_model"]
    assignments = pl.read_parquet(inp / cfg["fixed_assignments"]).filter(
        pl.col("strategy") == strategy
    )
    groups = (
        table.select("compound_id")
        .join(
            assignments.select("compound_id", "group"),
            on="compound_id",
            how="left",
            validate="1:1",
            maintain_order="left",
        )["group"]
        .fill_null("blind")
        .to_numpy()
    )
    scores, margins, predictions, decisions, memberships = [], [], [], [], []
    model_dir = out / "models"
    model_dir.mkdir()
    for seed in cfg["seeds"]:
        for bag, (train, holdout) in enumerate(
            folds(development, y, groups, strategy, cfg["bag_folds"], seed)
        ):
            assert set(groups[train]).isdisjoint(groups[holdout])
            assert not np.isin(train, blind).any()
            inner_score = np.full(len(table), np.nan)
            for inner, (fit, valid) in enumerate(
                folds(train, y, groups, strategy, cfg["inner_folds"], seed + bag + 1)
            ):
                assert set(groups[fit]).isdisjoint(groups[valid])
                assert not np.isin(np.r_[fit, valid], np.r_[holdout, blind]).any()
                fitted = fit_model(
                    kind,
                    fit,
                    fp,
                    x,
                    y,
                    tdi,
                    shift,
                    seed + bag * 10 + inner,
                    weights,
                    subtype,
                    tasks,
                    cfg,
                )
                inner_score[valid] = predict(kind, fitted, valid, fp, x)
            assert np.isfinite(inner_score[train]).all()
            cut = select_cut(y[train], inner_score[train], cfg["bootstrap_samples"])
            assert np.isfinite(cut)
            fitted = fit_model(
                kind, train, fp, x, y, tdi, shift, seed + bag * 10 + 9, weights, subtype, tasks, cfg
            )
            score = predict(kind, fitted, blind, fp, x)
            scores.append(score)
            margins.append(score - cut)
            model_path = model_dir / f"{seed}-bag{bag}.pt"
            torch.save(
                dict(
                    network_state=fitted["network"].state_dict(),
                    feature_mean=fitted["feature_mean"],
                    feature_scale=fitted["feature_scale"],
                    target_mean=fitted["means"],
                    target_scale=fitted["scales"],
                    config=cfg,
                    seed=seed,
                    fold=bag,
                    cut=cut,
                    architecture=[2048, 256, 64, tasks.shape[1]],
                ),
                model_path,
            )
            for row, s in zip(blind, score, strict=True):
                predictions.append(
                    dict(
                        compound_id=table["compound_id"][int(row)],
                        endpoint=cfg["endpoint"],
                        split_id=strategy,
                        fold=bag,
                        seed=seed,
                        model=kind,
                        score=float(s),
                        cut=cut,
                        margin=float(s - cut),
                    )
                )
            for row in train:
                memberships.append(
                    dict(
                        compound_id=table["compound_id"][int(row)],
                        endpoint=cfg["endpoint"],
                        split_id=strategy,
                        fold=bag,
                        seed=seed,
                        group=groups[row],
                        role="model_train",
                    )
                )
            decisions.append(
                dict(
                    seed=seed,
                    fold=bag,
                    cut=cut,
                    train_molecules=len(train),
                    train_classes=int(np.isfinite(y[train]).sum()),
                    model_path=str(model_path.name),
                )
            )
            print(json.dumps(decisions[-1]), flush=True)
    result = (
        table[blind]
        .select("compound_id", "molecule_name")
        .with_columns(
            pl.lit(cfg["endpoint"]).alias("endpoint"),
            pl.lit(strategy).alias("split_id"),
            pl.Series("score", np.mean(scores, axis=0)),
            pl.Series("margin", np.mean(margins, axis=0)),
            pl.Series("call", np.mean(margins, axis=0) > 0),
        )
    )
    assert result.height == 750 and not result["compound_id"].is_duplicated().any()
    result.write_parquet(out / "blind-policy.parquet")
    pl.DataFrame(predictions).write_parquet(out / "blind-fold-predictions.parquet")
    pl.DataFrame(memberships).write_parquet(out / "keyed-fit-memberships.parquet")
    pl.DataFrame(decisions).write_csv(out / "trained-cuts.csv")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                n_models=len(scores),
                positives=int(result["call"].sum()),
                models=decisions,
                no_upload=True,
                blind_assays_used=False,
                count_is_not_prescribed=True,
                limits="Same algorithm as nested validation; each deployed bag "
                "has 80 percent of development rather than 64 percent in outer validation. "
                "Generic upstream encoder corpus manifest is unknown. "
                "Weak-primary D6 class outcomes remain unavailable.",
            ),
            indent=2,
        )
        + "\n"
    )
