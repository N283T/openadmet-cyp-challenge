"""Refit the complete frozen September 15 recipe inside matched nested boundaries.

External encoder features are challenge-independent inputs. CYP graph weights,
readouts, PCA/scalers, graph epochs and decision cuts are trained afresh.
"""

import gc
import json
import os
import time
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop.data import BatchMolGraph
from evaluate_cyp_d6_zero_start import folds
from finetune_external_multitask import MultitaskModel, encoder_architecture
from lightgbm import LGBMClassifier
from rdkit import Chem
from scipy.special import expit
from scipy.stats import rankdata
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import ComplementNB
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from tabpfn import TabPFNClassifier
from tabpfn.constants import ModelVersion
from validate_cyp_d6_bag_policy import metric, select_cut


def compose_legacy(raw):
    """Exact raw-score core plus rank-transformed Morgan tails from the sent recipe."""
    assert raw.shape[1] == 11 and np.isfinite(raw).all()
    ranked = np.column_stack([rankdata(v, method="average") / (len(v) + 1) for v in raw.T])
    mixed_pfn = (raw[:, 5] + raw[:, 6]) / 2
    core_raw = np.mean(np.column_stack([raw[:, :5], mixed_pfn, raw[:, 7]]), axis=1)
    core = rankdata(core_raw, method="average") / (len(core_raw) + 1)
    live = (2 * core + ranked[:, 8]) / 3
    return ((2 * live + ranked[:, 9]) / 3 + (2 * live + ranked[:, 10]) / 3) / 2


class LegacyRecipe:
    def __init__(self, table, fp, tasks, cfg):
        self.table, self.fp, self.tasks, self.cfg = table, fp, tasks, cfg
        self.inp = Path(os.environ["CYP_INPUT_DIR"])
        ids = table.select("compound_id")

        def align(frame, prefix=None):
            cols = [
                c
                for c in frame.columns
                if c != "compound_id" and (prefix is None or c.startswith(prefix))
            ]
            values = ids.join(
                frame.select("compound_id", *cols),
                on="compound_id",
                how="left",
                validate="1:1",
                maintain_order="left",
            )
            result = values.select(cols).to_numpy().astype(np.float32)
            result[~np.isfinite(result)] = np.nan
            return result

        mordred = pl.concat(
            [
                pl.read_parquet(self.inp / f"mordred-{i:05d}.parquet")
                for i in range(cfg["mordred_files"])
            ]
        )
        xm = align(mordred)
        xc = align(pl.read_parquet(self.inp / "class-embeddings.parquet"), "emb_")
        mixed = align(pl.read_parquet(self.inp / "mixed-embeddings.parquet"), "emb_")
        assert xc.shape[1] == mixed.shape[1] == 600
        assert np.isfinite(xc).all() and np.isfinite(mixed).all()
        self.views = [fp, xm, np.column_stack([fp, xm]), xc, mixed]
        checkpoint = torch.load(
            self.inp / "external-best.pt", map_location="cpu", weights_only=False
        )
        self.encoder = {
            k.removeprefix("mp."): v
            for k, v in checkpoint["state_dict"].items()
            if k.startswith("mp.")
        }
        self.architecture = encoder_architecture(checkpoint, self.encoder)
        featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
        self.graphs = [featurizer(Chem.MolFromSmiles(s)) for s in table["standardized_smiles"]]
        self.fit_audit = []

    def graph_batch(self, indices):
        graph = BatchMolGraph([self.graphs[int(i)] for i in indices])
        graph.to("cuda")
        return graph

    def graph_predict(self, model, indices):
        model.eval()
        with torch.no_grad():
            return np.concatenate(
                [
                    model(self.graph_batch(indices[i : i + 128])).cpu().numpy()
                    for i in range(0, len(indices), 128)
                ]
            )

    @staticmethod
    def bce(logits, truth):
        mask = torch.isfinite(truth)
        error = torch.nn.functional.binary_cross_entropy_with_logits(
            logits, torch.nan_to_num(truth), reduction="none"
        )
        present = mask.any(0)
        return ((error * mask).sum(0) / mask.sum(0).clamp_min(1))[present].mean()

    def graph_fit(self, train, y, groups, strategy, seed):
        # Epoch selection sees only a molecular group holdout INSIDE this fit.
        fit, valid = folds(train, y, groups, strategy, 5, seed + 700)[0]
        assert set(groups[fit]).isdisjoint(groups[valid])
        relevant = np.isfinite(self.tasks[:, :2]).any(1)
        fit, valid = fit[relevant[fit]], valid[relevant[valid]]
        truth = self.tasks[:, :2].astype(np.float32)
        best_epoch = None
        for phase in ("epoch_selection", "refit"):
            rows = fit if phase == "epoch_selection" else train[relevant[train]]
            torch.manual_seed(seed)
            torch.cuda.manual_seed_all(seed)
            rng = np.random.default_rng(seed)
            model = MultitaskModel(self.encoder, 0.1, 2, self.architecture).to("cuda")
            optimizer = torch.optim.AdamW(
                [
                    {"params": model.mp.parameters(), "lr": 1e-5},
                    {"params": model.head.parameters(), "lr": 1e-4},
                ],
                weight_decay=1e-4,
            )
            best, stale = float("inf"), 0
            limit = self.cfg["graph_max_epochs"] if phase == "epoch_selection" else best_epoch
            for epoch in range(1, limit + 1):
                model.train()
                for batch in np.array_split(rng.permutation(rows), int(np.ceil(len(rows) / 64))):
                    loss = self.bce(
                        model(self.graph_batch(batch)), torch.tensor(truth[batch], device="cuda")
                    )
                    optimizer.zero_grad()
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(model.parameters(), 5)
                    optimizer.step()
                if phase == "epoch_selection":
                    value = float(
                        self.bce(
                            torch.tensor(self.graph_predict(model, valid)),
                            torch.tensor(truth[valid]),
                        )
                    )
                    if value < best - 1e-6:
                        best, best_epoch, stale = value, epoch, 0
                    else:
                        stale += 1
                    if stale >= self.cfg["graph_patience"]:
                        break
            if phase == "epoch_selection":
                del model, optimizer
        return model, int(best_epoch)

    def fit_predict(self, train, predict_sets, y, groups, strategy, seed, entry):
        start = time.monotonic()
        observed = train[np.isfinite(y[train])]
        assert len(np.unique(y[observed])) == 2
        for indices in predict_sets:
            assert not np.isin(train, indices).any()
            assert set(groups[train]).isdisjoint(groups[indices])
        scores = [[] for _ in predict_sets]
        for view in self.views:
            model = LGBMClassifier(
                objective="binary",
                n_estimators=400,
                num_leaves=15,
                min_child_samples=20,
                learning_rate=0.03,
                reg_lambda=2,
                n_jobs=3,
                random_state=seed,
                verbosity=-1,
                deterministic=True,
                force_col_wise=True,
            )
            model.fit(view[observed], y[observed].astype(int))
            for bank, indices in zip(scores, predict_sets, strict=True):
                bank.append(model.predict_proba(view[indices])[:, 1])
        mixed = self.views[-1]
        scaler = StandardScaler().fit(mixed[observed])
        scaled_train = scaler.transform(mixed[observed])
        pca = PCA(n_components=256, random_state=42).fit(scaled_train)
        for transform in (pca, None):
            train_x = transform.transform(scaled_train) if transform else scaled_train
            model = TabPFNClassifier.create_default_for_version(
                ModelVersion.V2_6,
                n_estimators=4,
                model_path=str(self.inp / "tabpfn-v2.6.ckpt"),
                device="cuda",
                random_state=seed,
                auto_scale_n_estimators=False,
                n_preprocessing_jobs=1,
                memory_saving_mode="auto",
                show_progress_bar=False,
                ignore_pretraining_limits=True,
                fit_mode="fit_preprocessors",
            )
            model.fit(train_x.astype(np.float32), y[observed].astype(int))
            positive = list(model.classes_).index(1)
            for bank, indices in zip(scores, predict_sets, strict=True):
                scaled = scaler.transform(mixed[indices])
                values = transform.transform(scaled) if transform else scaled
                bank.append(model.predict_proba(values.astype(np.float32))[:, positive])
            del model
            gc.collect()
            torch.cuda.empty_cache()
        model, epochs = self.graph_fit(train, y, groups, strategy, seed)
        for bank, indices in zip(scores, predict_sets, strict=True):
            bank.append(expit(self.graph_predict(model, indices)[:, 0]))
        del model
        for model in (
            LogisticRegression(C=0.05, max_iter=4000, solver="liblinear", random_state=seed),
            SVC(C=1.0, kernel="rbf"),
            ComplementNB(alpha=1.0),
        ):
            model.fit(self.fp[observed], y[observed].astype(int))
            for bank, indices in zip(scores, predict_sets, strict=True):
                bank.append(
                    model.decision_function(self.fp[indices])
                    if isinstance(model, SVC)
                    else model.predict_proba(self.fp[indices])[:, 1]
                )
        result = []
        for bank in scores:
            raw = np.column_stack(bank)
            result.append(compose_legacy(raw))
        self.fit_audit.append(
            {
                **entry,
                "fit_seed": seed,
                "fit_molecules": len(train),
                "fit_labels": len(observed),
                "graph_epochs": epochs,
                "seconds": time.monotonic() - start,
                "outer_test_used": False,
            }
        )
        print(json.dumps(self.fit_audit[-1]), flush=True)
        torch.cuda.empty_cache()
        return result


def evaluate(table, fp, x, y, tdi, shift, weights, subtype, tasks, development, blind, cfg):
    torch.set_num_threads(cfg["cpu_threads"])
    assert not np.isfinite(tasks[blind]).any()
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    recipe = LegacyRecipe(table, fp, tasks, cfg)
    assignments = pl.read_parquet(recipe.inp / cfg["fixed_assignments"])
    rows, cuts, audits, metrics = [], [], [], []
    for strategy in cfg["split_strategies"]:
        ordered = table.select("compound_id").join(
            assignments.filter(pl.col("strategy") == strategy).select(
                "compound_id", "fold", "group"
            ),
            on="compound_id",
            how="left",
            validate="1:1",
            maintain_order="left",
        )
        membership = ordered["fold"].to_numpy()
        groups = ordered["group"].fill_null("blind").to_numpy()
        assert np.isfinite(membership[development]).all()
        for seed in cfg["seeds"]:
            for outer in range(cfg["outer_folds"]):
                train = development[membership[development] != outer]
                test = development[membership[development] == outer]
                assert set(groups[train]).isdisjoint(groups[test])
                score_bank, margin_bank = [], []
                for bag, (bag_train, bag_valid) in enumerate(
                    folds(train, y, groups, strategy, cfg["bag_folds"], seed + outer * 100)
                ):
                    assert set(groups[bag_train]).isdisjoint(groups[bag_valid])
                    inner_score = np.full(len(table), np.nan)
                    for inner, (fit, valid) in enumerate(
                        folds(
                            bag_train,
                            y,
                            groups,
                            strategy,
                            cfg["inner_folds"],
                            seed + outer * 100 + bag + 1,
                        )
                    ):
                        assert not np.isin(np.r_[fit, valid], test).any()
                        inner_score[valid] = recipe.fit_predict(
                            fit,
                            [valid],
                            y,
                            groups,
                            strategy,
                            seed + outer * 1000 + bag * 10 + inner,
                            dict(
                                split_id=strategy,
                                seed=seed,
                                outer_fold=outer,
                                bag_fold=bag,
                                phase="cut_inner",
                                inner_fold=inner,
                            ),
                        )[0]
                    assert np.isfinite(inner_score[bag_train]).all()
                    cut = select_cut(y[bag_train], inner_score[bag_train], cfg["bootstrap_samples"])
                    score = recipe.fit_predict(
                        bag_train,
                        [test],
                        y,
                        groups,
                        strategy,
                        seed + outer * 1000 + bag * 10 + 9,
                        dict(
                            split_id=strategy,
                            seed=seed,
                            outer_fold=outer,
                            bag_fold=bag,
                            phase="bag_refit",
                            inner_fold=-1,
                        ),
                    )[0]
                    score_bank.append(score)
                    margin_bank.append(score - cut)
                    cuts.append(
                        dict(
                            split_id=strategy,
                            seed=seed,
                            outer_fold=outer,
                            model=cfg["model_name"],
                            bag_fold=bag,
                            threshold=cut,
                        )
                    )
                    pl.DataFrame(cuts).write_csv(out / "threshold-scale-audit.csv")
                    pl.DataFrame(recipe.fit_audit).write_csv(out / "fit-audit.csv")
                score, margin = np.mean(score_bank, axis=0), np.mean(margin_bank, axis=0)
                entry = dict(
                    split_id=strategy,
                    seed=seed,
                    fold=outer,
                    model=cfg["model_name"],
                    policy="centered_fold_mean",
                )
                metrics.append({**entry, **metric(y[test], score, margin > 0)})
                part = (
                    table[test]
                    .select("compound_id", "molecule_name")
                    .with_columns(
                        *[pl.lit(v).alias(k) for k, v in entry.items()],
                        pl.lit(cfg["endpoint"]).alias("endpoint"),
                        pl.Series("group", groups[test]),
                        pl.Series("label", y[test]),
                        pl.Series("score", score),
                        pl.Series("margin", margin),
                        pl.Series("call", margin > 0),
                    )
                )
                rows.append(part)
                pl.concat(rows).write_parquet(out / "nested-oof.parquet")
                pl.DataFrame(metrics).write_csv(out / "fold-metrics.csv")
                audits.append(
                    {
                        **entry,
                        "outer_test_used_in_fit_or_cut": False,
                        "mean_bag_threshold": float(np.mean([v["threshold"] for v in cuts[-5:]])),
                    }
                )
                pl.DataFrame(audits).write_csv(out / "leakage-and-operating-points.csv")
                print(json.dumps(metrics[-1]), flush=True)
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                metrics=metrics,
                fits=len(recipe.fit_audit),
                no_upload=True,
                limitations=(
                    "Frozen historical recipe evaluated prospectively on development; "
                    "selection of the historical recipe itself was not nested. Internal "
                    "graph epoch selection excludes each current fit holdout; no old CYP "
                    "weights/PCA/scalers reused. Score composition ranks each unlabelled "
                    "prediction batch, matching historical transductive rank composition. "
                    "Missing weak-primary D6 outcomes remain unidentified."
                ),
            ),
            indent=2,
        )
        + "\n"
    )
