"""Match the full-data submission timestamp and summarize immutable official CSVs."""

import json
import os
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
rows = []
for item in cfg["snapshots"]:
    table = pl.read_csv(inp / item["file"])
    ours = table.filter(pl.col("Username") == "N283T")
    baseline = table.filter(pl.col("Username") == "TabICL-baseline")
    assert ours.height == baseline.height == 1
    assert ours["Submitted"][0] == "2026-09-12 07:48 UTC"
    metric = "MA-ST-RAE" if item["endpoint"] == "MA" else "ST-RAE"
    rho = "MA-Spearman's ρ" if item["endpoint"] == "MA" else "Spearman's ρ"
    rows.append(
        dict(
            endpoint=item["endpoint"],
            rank=ours["Rank"][0],
            lb=ours[metric][0],
            official_tabicl=baseline[metric][0],
            difference=ours[metric][0] - baseline[metric][0],
            spearman=ours[rho][0],
            official_spearman=baseline[rho][0],
            oof=item["oof"],
            snapshot=item["snapshot"],
        )
    )
result = pl.DataFrame(rows)
result.write_csv(out / "comparison.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            submitted="2026-09-12 07:48 UTC",
            submission_id=cfg["submission_id"],
            timestamp_matched=True,
            scores=rows,
        ),
        indent=2,
    )
)
print(json.dumps(rows, indent=2))
