"""Package one legacy x gate D6 blind call vector as a classification candidate file.

The D6 score is rebuilt from the deployed Sep 15 blind score vector and the frozen
blind scores of run 7920b49a with the exact ranking of run 2d8d20f5 Part 3, and the
resulting top-k calls must reproduce that run's `blind_candidates.parquet` vector
for the configured variant. CYP3A4 is copied unchanged from the last sent payload
and must equal the historical best-scored Sep 15 vector. Nothing is fitted or
uploaded.
"""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.tracks import validate

IDENTIFIERS = ["SMILES", "Molecule_Name"]


def percentile_rank(values):
    """Identical to evaluate_cyp_d6_gate_concordance.percentile_rank."""
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values))
    ranks[order] = np.arange(len(values))
    return (ranks + 0.5) / len(values)


def top_calls(score, count):
    calls = np.zeros(len(score), dtype=bool)
    calls[np.argsort(-score, kind="stable")[:count]] = True
    return calls


def positions(score):
    return np.argsort(np.argsort(-score, kind="stable"), kind="stable") + 1


def variant_score(variant, legacy_rank, pool_rank, p_active):
    if variant == "L4_rankmean_legacy_pool":
        return 0.5 * (legacy_rank + pool_rank)
    if variant == "L2_legacy_x_sqrt_pa_morgan":
        return legacy_rank * np.sqrt(p_active)
    raise ValueError(f"Unsupported variant: {variant}")


def blind_column(test_scores, key):
    model, field = key.split(".")
    part = test_scores.filter(pl.col("model") == model).select(
        "compound_id", pl.col(field).alias(key)
    )
    assert part.height == 750 and part[key].null_count() == 0, key
    assert part[key].is_finite().all(), key
    return part


def payload(path, reference):
    frame = pl.read_parquet(path).sort("Molecule_Name")
    assert frame.height == 750
    assert frame.select(IDENTIFIERS).equals(reference.select(IDENTIFIERS))
    return frame


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    variant = cfg["variant"]
    spec = cfg["variants"][variant]
    d6, a4 = cfg["d6_endpoint"], cfg["a4_endpoint"]
    count = cfg["d6_positive_calls"]

    reference = pl.read_csv(inp / "reference.csv").select(IDENTIFIERS).sort("Molecule_Name")
    names = (
        pl.read_parquet(inp / "test_mapping.parquet")
        .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
        .unique()
    )
    assert names.height == 750 and names["compound_id"].n_unique() == 750
    last = payload(inp / "last_sent.parquet", reference)
    sep15 = payload(inp / "sep15_payload.parquet", reference)

    # A4: exact last-sent column, which must be the historical best-scored vector.
    assert last[a4].null_count() == 0 and last[a4].equals(sep15[a4])
    assert int(last[a4].sum()) == cfg["a4_positive_calls"]

    # Legacy Sep 15 blind score: must reproduce the 49 submitted Sep 15 D6 calls.
    deployed = pl.read_parquet(inp / "sep15_blind_scores.parquet").filter(
        pl.col("model") == cfg["legacy_model"]
    )
    assert deployed.height == 750 and deployed["score"].is_finite().all()
    table = deployed.select("compound_id", pl.col("score").alias("legacy_score")).join(
        names, on="compound_id", how="left", validate="1:1"
    )
    test_scores = pl.read_parquet(inp / "test_scores.parquet")
    for key in (cfg["pool_column"], cfg["gate_column"]):
        table = table.join(blind_column(test_scores, key), on="compound_id", validate="1:1")
    # Same compound order as run 2d8d20f5 Part 3 so stable tie handling matches.
    table = table.sort("compound_id")
    assert table.height == 750 and table["Molecule_Name"].null_count() == 0
    status = (
        table.select("Molecule_Name")
        .join(
            sep15.select("Molecule_Name", pl.col(d6).alias("sep15_call")),
            on="Molecule_Name",
            validate="1:1",
            maintain_order="left",
        )
        .join(
            last.select("Molecule_Name", pl.col(d6).alias("last_sent_call")),
            on="Molecule_Name",
            validate="1:1",
            maintain_order="left",
        )
    )
    sep15_call = status["sep15_call"].to_numpy()
    last_call = status["last_sent_call"].to_numpy()
    legacy = table["legacy_score"].to_numpy()
    assert np.array_equal(top_calls(legacy, int(sep15_call.sum())), sep15_call)
    assert int(sep15_call.sum()) == count

    legacy_rank = percentile_rank(legacy)
    pool_rank = percentile_rank(table[cfg["pool_column"]].to_numpy())
    p_active = table[cfg["gate_column"]].to_numpy()
    score = variant_score(variant, legacy_rank, pool_rank, p_active)
    calls = top_calls(score, count)
    order = np.argsort(-score, kind="stable")
    boundary_gap = float(score[order[count - 1]] - score[order[count]])
    assert boundary_gap > 0, "Tied call boundary"

    # Cross-check against the reviewed blind vector of the selection run.
    reviewed = (
        pl.read_parquet(inp / "blind_candidates.parquet")
        .filter((pl.col("variant") == variant) & (pl.col("n_calls") == count))
        .rename({"molecule_name": "Molecule_Name"})
    )
    assert reviewed.height == 750
    reviewed = table.select("Molecule_Name").join(
        reviewed.select("Molecule_Name", "call"),
        on="Molecule_Name",
        validate="1:1",
        maintain_order="left",
    )
    assert np.array_equal(reviewed["call"].to_numpy(), calls)

    swaps = dict(
        shared_sep15=int((calls & sep15_call).sum()),
        added=int((calls & ~sep15_call).sum()),
        removed=int((~calls & sep15_call).sum()),
    )
    assert swaps == spec["expected"], (swaps, spec["expected"])
    versus_last = dict(
        last_sent_positive_calls=int(last_call.sum()),
        shared=int((calls & last_call).sum()),
        added=int((calls & ~last_call).sum()),
        removed=int((~calls & last_call).sum()),
        changed=int((calls != last_call).sum()),
    )

    d6_calls = table.select("Molecule_Name").with_columns(pl.Series(d6, calls))
    candidate = (
        last.drop(d6)
        .join(d6_calls, on="Molecule_Name", validate="1:1", maintain_order="left")
        .select(*IDENTIFIERS, d6, a4)
    )
    assert candidate[a4].equals(last[a4]) and candidate[a4].equals(sep15[a4])
    assert int(candidate[d6].sum()) == count
    candidate.write_parquet(out / "predictions.parquet")
    candidate.write_csv(out / "predictions.csv")
    validation = validate(out / "predictions.parquet", "classification", inp / "reference.csv")

    scores = table.select("compound_id", "Molecule_Name").with_columns(
        pl.lit(d6).alias("endpoint"),
        pl.lit("blind").alias("split_id"),
        pl.lit(-1).alias("fold"),
        pl.lit(variant).alias("model"),
        pl.Series("legacy_score", legacy),
        pl.Series("legacy_rank", legacy_rank),
        pl.Series("pool_score", table[cfg["pool_column"]].to_numpy()),
        pl.Series("pool_rank", pool_rank),
        pl.Series("p_active_morgan", p_active),
        pl.Series("score", score),
        pl.Series("position", positions(score)),
        pl.Series("legacy_position", positions(legacy)),
        pl.Series("call", calls),
        pl.Series("sep15_call", sep15_call),
        pl.Series("last_sent_call", last_call),
    )
    scores.write_parquet(out / "test_scores.parquet")
    changes = (
        scores.filter(pl.col("call") != pl.col("sep15_call"))
        .with_columns(
            pl.when(pl.col("call"))
            .then(pl.lit("added"))
            .otherwise(pl.lit("removed"))
            .alias("change")
        )
        .sort("change", "position")
    )
    changes.write_csv(out / "swaps_vs_sep15.csv")

    summary = dict(
        variant=variant,
        role=spec["role"],
        formula=spec["formula"],
        selection_run=cfg["selection_run"],
        score_run=cfg["score_run"],
        legacy_run=cfg["legacy_run"],
        sep15_submission=cfg["sep15_submission"],
        last_sent_submission=cfg["last_sent_submission"],
        d6_positive_calls=int(calls.sum()),
        d6_boundary_gap=boundary_gap,
        d6_vs_sep15=swaps,
        d6_vs_last_sent=versus_last,
        d6_matches_selection_run_vector=True,
        mean_p_active_called=float(p_active[calls].mean()),
        mean_p_active_sep15=float(p_active[sep15_call].mean()),
        a4_positive_calls=int(candidate[a4].sum()),
        a4_changed_vs_last_sent=0,
        a4_changed_vs_sep15=0,
        validation=validation,
        uploaded=False,
        prediction_sha256=hashlib.sha256((out / "predictions.parquet").read_bytes()).hexdigest(),
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    lines = [
        f"# D6 gate candidate package: `{variant}` ({spec['role']})",
        "",
        f"D6 score: `{spec['formula']}` over the 750 blind compounds; top {count} called.",
        f"Legacy = deployed Sep 15 blind vector (`{cfg['legacy_model']}`, run "
        f"`{cfg['legacy_run']}`); pool and gate from run `{cfg['score_run']}`.",
        f"The call vector reproduces run `{cfg['selection_run']}` `blind_candidates.parquet`.",
        "",
        f"D6 vs Sep 15 payload: {swaps['shared_sep15']} shared, {swaps['added']} added, "
        f"{swaps['removed']} removed (boundary gap {boundary_gap:.6f}).",
        f"D6 vs last sent payload `{cfg['last_sent_submission']}` "
        f"({versus_last['last_sent_positive_calls']} calls): {versus_last['shared']} shared, "
        f"{versus_last['added']} added, {versus_last['removed']} removed.",
        f"Mean Morgan gate p_active of called compounds {summary['mean_p_active_called']:.3f} "
        f"(Sep 15: {summary['mean_p_active_sep15']:.3f}).",
        f"A4: {summary['a4_positive_calls']} calls, bitwise identical to the last sent payload "
        "and to the historical best-scored Sep 15 vector.",
        "",
        "## D6 changes versus the Sep 15 payload",
        "",
        "| Change | Molecule | Position | Legacy position | Pool rank | p_active |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in changes.to_dicts():
        lines.append(
            f"| {row['change']} | {row['Molecule_Name']} | {row['position']} | "
            f"{row['legacy_position']} | {row['pool_rank']:.3f} | {row['p_active_morgan']:.3f} |"
        )
    lines += ["", "No model was fitted and nothing was uploaded.", ""]
    (out / "report.md").write_text("\n".join(lines))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
