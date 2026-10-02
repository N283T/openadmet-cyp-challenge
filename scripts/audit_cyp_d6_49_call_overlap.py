"""Compare exact submitted model top-49 identities with the historical best payload."""

import json
import os
from pathlib import Path

import polars as pl


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    endpoint, k = cfg["endpoint"], cfg["diagnostic_top_k"]
    policy = pl.read_parquet(inp / "blind-policy.parquet").rename(
        {"molecule_name": "Molecule_Name"}
    )
    best = pl.read_parquet(inp / "best.parquet").sort("Molecule_Name")
    sent = pl.read_parquet(inp / "submitted.parquet").sort("Molecule_Name")
    assert best.select("Molecule_Name", "SMILES").equals(sent.select("Molecule_Name", "SMILES"))
    assert policy.height == 750 and policy["Molecule_Name"].n_unique() == 750
    assert best[endpoint].sum() == k and sent[endpoint].sum() == cfg["submitted_positives"]
    ranked = policy.sort(["margin", "Molecule_Name"], descending=[True, False]).with_row_index(
        "new_rank", offset=1
    )
    ranked = ranked.join(
        best.select("Molecule_Name", pl.col(endpoint).alias("best_call")),
        on="Molecule_Name",
        validate="1:1",
        maintain_order="left",
    ).join(
        sent.select("Molecule_Name", pl.col(endpoint).alias("sent_call")),
        on="Molecule_Name",
        validate="1:1",
        maintain_order="left",
    )
    assert not ranked["best_call"].has_nulls() and not ranked["sent_call"].has_nulls()
    assert ranked["call"].equals(ranked["sent_call"])
    assert ((ranked["margin"] > 0) == ranked["sent_call"]).all()
    assert set(ranked.head(k)["Molecule_Name"]) == set(
        policy.sort(["score", "Molecule_Name"], descending=[True, False]).head(k)["Molecule_Name"]
    )
    ranked = ranked.with_columns((pl.col("new_rank") <= k).alias("top49_call"))
    rows = []
    for name, column in [("counterfactual_top49", "top49_call"), ("submitted175", "sent_call")]:
        rows.append(
            dict(
                policy=name,
                positives=int(ranked[column].sum()),
                shared=int(ranked.filter(pl.col(column) & pl.col("best_call")).height),
                added=int(ranked.filter(pl.col(column) & ~pl.col("best_call")).height),
                removed=int(ranked.filter(~pl.col(column) & pl.col("best_call")).height),
                changed=int(ranked.filter(pl.col(column) != pl.col("best_call")).height),
            )
        )
    ranked.write_parquet(out / "keyed-call-comparison.parquet")
    ranked.filter(pl.col("best_call")).write_csv(out / "historical-best-new-ranks.csv")
    ranked.filter(pl.col("top49_call") != pl.col("best_call")).write_csv(out / "top49-swaps.csv")
    pl.DataFrame(rows).write_csv(out / "overlap-summary.csv")
    (out / "summary.json").write_text(
        json.dumps(
            dict(config=cfg, comparisons=rows, no_upload=True, no_lb_estimate=True), indent=2
        )
        + "\n"
    )
    print(json.dumps(rows))


if __name__ == "__main__":
    main()
