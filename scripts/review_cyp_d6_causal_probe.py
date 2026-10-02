"""Connect historical failures to decision evidence and freeze an independent probe."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, matthews_corrcoef

from cyp_challenge.tracks import validate


def calls_at_rate(score, folds, rate):
    calls = np.zeros(len(score), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        order = np.argsort(-score[rows], kind="stable")
        calls[rows[order[: int(round(rate * len(rows)))]]] = True
    return calls


def projected(y, calls, prior):
    tpr, fpr = calls[y == 1].mean(), calls[y == 0].mean()
    tp, fn = prior * tpr, prior * (1 - tpr)
    fp, tn = (1 - prior) * fpr, (1 - prior) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    endpoint = cfg["endpoint"]
    new = pl.read_parquet(inp / "new-oof.parquet")
    threshold = pl.read_csv(inp / "new-operating-points.csv")
    summary = json.loads((inp / "new-summary.json").read_text())
    all_curves, bank = [], []
    for split in ["umap", "scaffold"]:
        old = pl.read_parquet(inp / f"control-{split}.parquet")
        ap = pl.read_parquet(inp / f"ap-{split}.parquet")
        historical = pl.concat([old, ap]).select(
            "compound_id",
            "fold",
            "group",
            "split_id",
            "model",
            "score",
            pl.col("y_true").cast(pl.Float64).alias("label"),
        )
        modern = new.filter(pl.col("split_id") == split).select(
            "compound_id", "fold", "split_id", "model", "score", "label"
        )
        modern = modern.join(
            old.filter(pl.col("model") == "sep15_control_raw_oof_analog").select(
                "compound_id", "group"
            ),
            on="compound_id",
            how="left",
        )
        data = pl.concat([historical, modern.select(historical.columns)]).sort(
            "model", "compound_id"
        )
        bank.append(data)
        for scope in ["full_pool", "labelled_only"]:
            scoped = data if scope == "full_pool" else data.filter(pl.col("label").is_finite())
            reference = scoped.filter(pl.col("model") == "sep15_control_raw_oof_analog").sort(
                "compound_id"
            )
            for model in scoped["model"].unique().sort():
                frame = scoped.filter(pl.col("model") == model).sort("compound_id")
                assert frame["compound_id"].equals(reference["compound_id"])
                assert frame["fold"].equals(reference["fold"])
                np.testing.assert_allclose(
                    frame["label"].to_numpy(), reference["label"].to_numpy(), equal_nan=True
                )
                labels, score, folds = (
                    frame["label"].to_numpy(),
                    frame["score"].to_numpy(),
                    frame["fold"].to_numpy(),
                )
                scored = np.isfinite(labels)
                y = labels[scored].astype(int)
                for rate in cfg["rates"]:
                    calls = calls_at_rate(score, folds, rate)
                    control = calls_at_rate(reference["score"].to_numpy(), folds, rate)
                    added, removed = calls & ~control, ~calls & control
                    for prior in cfg["prior_sensitivity"]:
                        all_curves.append(
                            dict(
                                split_id=split,
                                scope=scope,
                                model=model,
                                rate=rate,
                                prior=prior,
                                pool_rows=len(frame),
                                scored_rows=int(scored.sum()),
                                calls_all=int(calls.sum()),
                                calls_scored=int(calls[scored].sum()),
                                mcc=float(matthews_corrcoef(y, calls[scored])),
                                projected_mcc=projected(y, calls[scored], prior),
                                delta_control=projected(y, calls[scored], prior)
                                - projected(y, control[scored], prior),
                                added_labelled=int((added & scored).sum()),
                                removed_labelled=int((removed & scored).sum()),
                                added_positive=int((added & (labels == 1)).sum()),
                                removed_positive=int((removed & (labels == 1)).sum()),
                            )
                        )
    pl.concat(bank).write_parquet(out / "matched-model-split-bank.parquet")
    curves = pl.DataFrame(all_curves)
    curves.write_csv(out / "full-pool-decision-curves.csv")
    # Retrospective October 1 objective mismatch: scores and decisions are fixed,
    # and were conditional on the labelled pool, not complete-pool predictions.
    trunk = pl.read_parquet(inp / "trunk-oof.parquet").filter(pl.col("endpoint") == endpoint)
    trunk_rows = []
    for model in ["current", "single_selected_AP", "multi_block_selected_AP"]:
        frame = trunk.filter(pl.col("model") == model).sort("compound_id")
        y, calls = frame["y_true"].to_numpy(), frame["call"].to_numpy()
        trunk_rows.append(
            dict(
                model=model,
                n=len(frame),
                ap=float(average_precision_score(y, frame["prediction"].to_numpy())),
                ordinary_mcc=float(matthews_corrcoef(y, calls)),
                projected_mcc_007=projected(y, calls, 0.07),
            )
        )
    pl.DataFrame(trunk_rows).write_csv(out / "failed-trunk-objective-audit.csv")
    fold_scores = pl.read_parquet(inp / "new-blind-folds.parquet").join(
        threshold.filter(pl.col("fold") >= 0).select("split_id", "fold", "model", "threshold"),
        on=["split_id", "fold", "model"],
        how="left",
        validate="m:1",
    )
    assert fold_scores["threshold"].is_null().sum() == 0
    inference = []
    for split in ["umap", "scaffold"]:
        for model in summary["config"]["models"]:
            fold = fold_scores.filter((pl.col("split_id") == split) & (pl.col("model") == model))
            reduced = fold.group_by("compound_id").agg(
                pl.col("score").mean().alias("fold_mean"),
                (pl.col("score") - pl.col("threshold")).mean().alias("mean_margin"),
                (pl.col("score") > pl.col("threshold")).sum().alias("positive_votes"),
            )
            full = pl.read_parquet(inp / f"blind-{split}-{model}.parquet")
            joined = reduced.join(
                full.select(
                    "compound_id",
                    pl.col("score").alias("full_score"),
                    pl.col("call").alias("full_call"),
                ),
                on="compound_id",
                validate="1:1",
            )
            inference.append(
                dict(
                    split_id=split,
                    model=model,
                    full_refit_calls=int(joined["full_call"].sum()),
                    centered_fold_calls=int((joined["mean_margin"] > 0).sum()),
                    majority_vote_calls=int((joined["positive_votes"] >= 3).sum()),
                    full_vs_centered_changed=int(
                        (joined["full_call"] != (joined["mean_margin"] > 0)).sum()
                    ),
                    full_vs_fold_rank_correlation=float(
                        np.corrcoef(rankdata(joined["full_score"]), rankdata(joined["fold_mean"]))[
                            0, 1
                        ]
                    ),
                )
            )
            joined.write_parquet(out / f"inference-{split}-{model}.parquet")
    pl.DataFrame(inference).write_csv(out / "inference-sensitivity.csv")
    chosen = fold_scores.filter(
        (pl.col("model") == cfg["candidate_model"]) & (pl.col("split_id") == cfg["candidate_split"])
    )
    candidate_scores = chosen.group_by("compound_id").agg(
        (pl.col("score") - pl.col("threshold")).mean().alias("margin"),
        (pl.col("score") > pl.col("threshold")).sum().alias("positive_votes"),
    )
    assert len(candidate_scores) == 750 and chosen.group_by("compound_id").len()["len"].min() == 5
    aliases = (
        pl.read_parquet(inp / "aliases.parquet")
        .filter(pl.col("role") == "test")
        .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
        .unique()
    )
    last = pl.read_parquet(inp / "last-sent.parquet").sort("Molecule_Name")
    best = pl.read_parquet(inp / "best-control.parquet").sort("Molecule_Name")
    keyed = (
        last.join(aliases, on="Molecule_Name", validate="1:1")
        .join(candidate_scores, on="compound_id", validate="1:1")
        .sort("Molecule_Name")
    )
    candidate = keyed.select(last.columns).with_columns(pl.Series(endpoint, keyed["margin"] > 0))
    candidate.write_parquet(out / "predictions.parquet")
    candidate.write_csv(out / "predictions.csv")
    keyed.with_columns(pl.Series("candidate_call", keyed["margin"] > 0)).write_parquet(
        out / "blind-review.parquet"
    )
    validation = validate(out / "predictions.parquet", "classification", inp / "reference.csv")
    assert candidate["CYP3A4_is_TDI"].equals(last["CYP3A4_is_TDI"])
    changes = []
    for name, control in [("last_sent", last), ("best_D6", best)]:
        for target in [endpoint, "CYP3A4_is_TDI"]:
            old, new_call = control[target].to_numpy(), candidate[target].to_numpy()
            changes.append(
                dict(
                    reference=name,
                    endpoint=target,
                    positives=int(new_call.sum()),
                    added=int((new_call & ~old).sum()),
                    removed=int((~new_call & old).sum()),
                    retained=int((new_call & old).sum()),
                )
            )
    pl.DataFrame(changes).write_csv(out / "target-changes.csv")
    live = json.loads((inp / "live-summary.json").read_text())
    result = dict(
        config=cfg,
        validation=validation,
        training_summary=summary,
        failed_trunk_objective=trunk_rows,
        inference=inference,
        target_changes=changes,
        hypotheses=[
            "October 1 D6 loss changes called identities at fixed count, not only rate",
            "AP gain did not imply decision gain in labelled-only trunk diagnostics",
            "Earlier complete-pool UMAP and scaffold still missed the AP-mass90 Live reversal",
            "New network has inference/threshold sensitivity; this does not explain prior losses",
        ],
        prospective_scope="Independent-model/variable-decision probe; no causal isolation",
        limitations=[
            "Unknown individual blind labels and scored mask",
            "No proven Live gain",
            "TPR/FPR prior projections assume conditional transport that may fail",
            "Blind five-fold averaging is not directly validated by single-model outer OOF calls",
            "Candidate model/policy chosen after development inspection",
        ],
        historical_live_summary=live,
    )
    (out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    report = [
        "# D6 failure-linked comparison and prospective probe",
        "",
        "The candidate tests independent D6 learning, inner-selected probability operating "
        "points and centered fold inference. It uses no historical score mixture and no prescribed "
        "positive count. A4 is copied from the last sent payload.",
        "",
        "## Matched historical splits",
        "",
        "| Model | Split | AP | Nested bootstrap MCC |",
        "| --- | --- | ---: | ---: |",
    ]
    for row in summary["aggregate"]:
        report.append(
            f"| {row['model']} | {row['split_id']} | {row['ap']:.4f} | {row['bootstrap_mcc']:.4f} |"
        )
    report += [
        "",
        "## Connection to prior failures",
        "",
        "October 1 held the D6 call count at 49 while exchanging 21 calls each way. Live MCC "
        "fell .308635 to .104986. The called identities changed at fixed count.",
        "",
        "The conditional multi-block AP strategy linked to the probe improved AP, "
        "but its historical-rate projected MCC "
        "did not improve. The AP-mass90 counterexample survived full-pool and scaffold "
        "checks, so neither correcting pool scope nor changing splits alone explains all failures. "
        "Full-pool and labelled-only curves are retained as sensitivity evidence; "
        "missing labels stay missing. Historical 6.533% is only a retrospective comparison.",
        "The conditional strategy result is not an independent validation of the exact "
        "full-development-selected blind payload. The single-model strategy is shown separately.",
        "",
        "## Inference sensitivity",
        "",
        "| Model | Split | Full refit calls | Centered-fold calls | Changed |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    for row in inference:
        report.append(
            f"| {row['model']} | {row['split_id']} | {row['full_refit_calls']} | "
            f"{row['centered_fold_calls']} | {row['full_vs_centered_changed']} |"
        )
    report += ["", "## Reviewed payload", ""]
    for row in changes:
        report.append(
            f"- Against {row['reference']}, {row['endpoint']}: {row['positives']} positives, "
            f"{row['added']} additions, {row['removed']} removals, {row['retained']} retained."
        )
    report += [
        "",
        "## Interpretation and limits",
        "",
        "This is a submission candidate for a prospective experiment. It does not establish "
        "that the historic reversals are solved. The five-member blind aggregation is not an "
        "independent nested-CV evaluation of an identical bag; it is an inference choice whose "
        "sensitivity is explicitly reported. One Live observation tests the combined new pipeline "
        "and cannot isolate one cause. No individual hidden label is reconstructed. "
        "Uncertainty is documented for the submission decision rather than an automatic veto.",
        "",
    ]
    (out / "report.md").write_text("\n".join(report))
    print(
        json.dumps(
            dict(target_changes=changes, failed_trunk_objective=trunk_rows, inference=inference),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
