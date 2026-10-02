"""Inspect and freeze a D6 AP rank-blend probe without uploading."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata, spearmanr

from cyp_challenge.tracks import validate

D6, A4 = "CYP2D6_is_TDI", "CYP3A4_is_TDI"


def ranks(score):
    return rankdata(score, method="average") / (len(score) + 1)


def top(score, count):
    calls = np.zeros(len(score), dtype=bool)
    calls[np.argsort(-score, kind="stable")[:count]] = True
    return calls


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    decision = json.loads((inp / "decision.json").read_text())
    if decision["kind"] != "blend":
        raise ValueError("This review packager requires the selected rank-blend family")
    base = pl.read_parquet(inp / "last-sent.parquet").sort("Molecule_Name")
    best = pl.read_parquet(inp / "performance-control.parquet").sort("Molecule_Name")
    reference = pl.read_csv(inp / "reference.csv").sort("Molecule_Name")
    aliases = (
        pl.read_parquet(inp / "aliases.parquet")
        .filter(pl.col("role") == "test")
        .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
        .unique()
    )
    ids = reference.join(aliases, on="Molecule_Name", validate="1:1", maintain_order="left")
    assert ids.height == 750 and ids["compound_id"].n_unique() == 750
    assert ids["compound_id"].null_count() == 0
    old_scores = pl.read_parquet(inp / "control-scores.parquet")
    ids = ids.join(
        old_scores.select("compound_id", "d6_score", "d6_call"),
        on="compound_id",
        validate="1:1",
        maintain_order="left",
    )
    assert ids["d6_score"].is_finite().all()
    control_rank = ranks(ids["d6_score"].to_numpy())
    count = decision["positive_calls"]
    old_call = top(control_rank, count)
    assert np.array_equal(old_call, ids["d6_call"].to_numpy())
    assert np.array_equal(old_call, base[D6].to_numpy())
    assert base.select("SMILES", "Molecule_Name").equals(
        reference.select("SMILES", "Molecule_Name")
    )
    assert best.select("SMILES", "Molecule_Name").equals(
        reference.select("SMILES", "Molecule_Name")
    )
    assert best[D6].equals(base[D6])
    raw = pl.read_parquet(inp / "trunk-predictions.parquet").filter(
        (pl.col("role") == "test") & (pl.col("endpoint") == D6)
    )
    fold_blocks = []
    mean_blocks = []
    prediction_frames = []
    for model in decision["models"]:
        part = raw.filter(pl.col("model") == model)
        assert part.height == 3750 and sorted(part["fold"].unique()) == list(range(5))
        fold_scores = []
        for fold in range(5):
            aligned = ids.select("compound_id").join(
                part.filter(pl.col("fold") == fold).select("compound_id", "prediction"),
                on="compound_id",
                validate="1:1",
                maintain_order="left",
            )
            assert aligned["prediction"].is_finite().all()
            values = aligned["prediction"].to_numpy()
            fold_scores.append(values)
            prediction_frames.append(
                ids.select("compound_id").with_columns(
                    pl.lit(D6).alias("endpoint"),
                    pl.lit("blind").alias("split_id"),
                    pl.lit(fold).alias("fold"),
                    pl.lit(model).alias("model"),
                    pl.Series("prediction", values),
                )
            )
        means = np.column_stack(fold_scores).mean(axis=1)
        mean_blocks.append(ranks(means))
        fold_blocks.append(np.column_stack([ranks(values) for values in fold_scores]))
    trunk_rank = np.mean(mean_blocks, axis=0)
    weight = decision["weight"]
    blend = (1 - weight) * control_rank + weight * trunk_rank
    calls = top(blend, count)
    candidate = base.with_columns(pl.Series(D6, calls))
    assert candidate[A4].equals(base[A4])
    candidate.write_parquet(out / "predictions.parquet")
    candidate.write_csv(out / "predictions.csv")
    validation = validate(out / "predictions.parquet", "classification", inp / "reference.csv")
    votes = np.zeros(750, dtype=int)
    fold_call_rows = []
    for fold in range(5):
        fold_trunk_rank = np.mean([block[:, fold] for block in fold_blocks], axis=0)
        fold_blend = (1 - weight) * control_rank + weight * fold_trunk_rank
        fold_calls = top(fold_blend, count)
        votes += fold_calls
        fold_call_rows.append(
            {
                "fold": fold,
                "positive_calls": int(fold_calls.sum()),
                "differences_from_mean_candidate": int((fold_calls != calls).sum()),
            }
        )
    neighbor = pl.read_parquet(inp / "neighbor-diagnostics.parquet").select(
        "compound_id",
        "nearest_any_similarity",
        "nearest_positive_similarity",
        "nearest_negative_similarity",
        "top10_similarity_weighted_positive_rate",
    )
    score_table = (
        ids.select("compound_id", "Molecule_Name", "SMILES")
        .with_columns(
            pl.lit(D6).alias("endpoint"),
            pl.lit("blind").alias("split_id"),
            pl.lit(-1).alias("fold"),
            pl.lit("full_development_AP_rank_blend").alias("model"),
            pl.Series("control_rank", control_rank),
            pl.Series("trunk_rank", trunk_rank),
            pl.Series("prediction", blend),
            pl.Series("control_call", old_call),
            pl.Series("call", calls),
            pl.Series(
                "control_position",
                np.argsort(np.argsort(-control_rank, kind="stable"), kind="stable") + 1,
            ),
            pl.Series(
                "candidate_position",
                np.argsort(np.argsort(-blend, kind="stable"), kind="stable") + 1,
            ),
            pl.Series("trunk_fold_positive_votes", votes),
        )
        .join(neighbor, on="compound_id", validate="1:1", maintain_order="left")
    )
    assert score_table["nearest_any_similarity"].is_finite().all()
    changes = (
        score_table.filter(pl.col("call") != pl.col("control_call"))
        .with_columns(
            pl.when(pl.col("call"))
            .then(pl.lit("added"))
            .otherwise(pl.lit("removed"))
            .alias("change")
        )
        .sort("change", "candidate_position")
    )
    score_table.write_parquet(out / "test-scores.parquet")
    changes.write_csv(out / "blind-swaps.csv")
    pl.concat(prediction_frames).write_parquet(out / "selected-trunk-fold-test.parquet")
    neighbors = (
        changes.group_by("change")
        .agg(
            pl.len().alias("n"),
            pl.col("top10_similarity_weighted_positive_rate")
            .mean()
            .alias("mean_neighbor_positive_rate"),
            pl.col("nearest_positive_similarity").mean().alias("mean_nearest_positive_similarity"),
            pl.col("nearest_negative_similarity").mean().alias("mean_nearest_negative_similarity"),
            pl.col("trunk_fold_positive_votes").mean().alias("mean_trunk_fold_positive_votes"),
        )
        .sort("change")
    )
    neighbors.write_csv(out / "swap-summary.csv")
    sensitivity = []
    for w in [0.0, 0.1, 0.25, weight, 1.0]:
        c = top((1 - w) * control_rank + w * trunk_rank, count)
        sensitivity.append(
            {
                "trunk_weight": w,
                "positive_calls": int(c.sum()),
                "added": int((c & ~old_call).sum()),
                "removed": int((~c & old_call).sum()),
                "retained": int((c & old_call).sum()),
            }
        )
    pl.DataFrame(sensitivity).unique(maintain_order=True).write_csv(out / "weight-sensitivity.csv")
    pl.DataFrame(fold_call_rows).write_csv(out / "trunk-fold-call-sensitivity.csv")
    swaps = {
        "added": int((calls & ~old_call).sum()),
        "removed": int((~calls & old_call).sum()),
        "retained": int((calls & old_call).sum()),
    }
    summary = {
        "selection_run": cfg["selection_run"],
        "decision": decision,
        "last_sent_submission": cfg["last_sent_submission"],
        "performance_control_submission": cfg["performance_control_submission"],
        "d6_positive_calls": int(calls.sum()),
        "d6_swaps": swaps,
        "a4_positive_calls": int(candidate[A4].sum()),
        "a4_changed_vs_last_sent": 0,
        "a4_changed_vs_performance_control": int((candidate[A4] != best[A4]).sum()),
        "score_spearman_vs_control": float(spearmanr(blend, control_rank).statistic),
        "trunk_score_spearman_vs_control": float(spearmanr(trunk_rank, control_rank).statistic),
        "neighbor_summary": neighbors.to_dicts(),
        "validation": validation,
        "uploaded": False,
        "prediction_sha256": hashlib.sha256((out / "predictions.parquet").read_bytes()).hexdigest(),
        "limitations": (
            "AP-ranked exploratory probe; incomplete-pool conditional OOF. "
            "Neighbor summaries are descriptive, not hidden-label estimates. "
            "Fold sensitivity varies the trunk inference only "
            "with the original blind control fixed."
        ),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = [
        "# Unsubmitted D6 AP-trunk probe: blind call review",
        "",
        "The user authorized inspection only and explicitly prohibited submission.",
        "",
        f"Selected full-development recipe: `{decision['models']}`; "
        f"trunk rank weight {weight:.2f}.",
        "Core-pocket, Fe and remaining HEM pair embeddings "
        "are aggregated separately and read with LightGBM. "
        "Five raw fold predictions are averaged, then ranked across the 750 blind compounds. "
        "The final score blends this rank with the existing D6 ensemble rank.",
        "",
        f"D6: {count} positives; {swaps['retained']} retained, "
        f"{swaps['added']} added, {swaps['removed']} removed.",
        f"A4: {summary['a4_positive_calls']} positives, "
        "unchanged from the last sent 275-call probe. "
        "This does not restore the earlier best-performing 249-call vector.",
        "",
        "The previously reported AP .296860→.318808 is the conditional held-fold procedure result, "
        "not a validation score measured for this exact final blind vector. The full-development "
        "inner AP is a recipe-selection diagnostic, not independent confirmation.",
        "",
        "## Changed blind calls",
        "",
        "| Change | Molecule | Control position | Candidate position | "
        "Weighted positive-neighbor fraction | Trunk-fold positive votes |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in changes.to_dicts():
        lines.append(
            f"| {row['change']} | {row['Molecule_Name']} | {row['control_position']} | "
            f"{row['candidate_position']} | {row['top10_similarity_weighted_positive_rate']:.3f} | "
            f"{row['trunk_fold_positive_votes']}/5 |"
        )
    lines.extend(
        [
            "",
            "Full keyed scores, SMILES, nearest positive/negative similarities, weaker-weight "
            "sensitivity and trunk-fold sensitivity are retained alongside this report.",
            "",
            "No upload or model training was performed. No full-pool features were added.",
            "",
        ]
    )
    (out / "report.md").write_text("\n".join(lines))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
