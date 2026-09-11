"""Register audited series-CV metrics with immutable input provenance."""

import json
import os
from pathlib import Path

import polars as pl
import psycopg

from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
count = 0
with psycopg.connect(**cfg["database"]) as conn:
    for i, parent in enumerate(cfg["parents"]):
        assert (
            conn.execute("SELECT status FROM runs WHERE id=%s", (parent,)).fetchone()[0]
            == "completed"
        )
        f = pl.read_csv(next(inp.glob(f"{i:03d}-*")))
        for r in f.iter_rows(named=True):
            for metric in ["ST_RAE", "MAE", "Spearman"]:
                conn.execute(
                    "INSERT INTO metrics VALUES (%s,%s,%s,%s,%s,%s)",
                    (
                        Path(os.environ["CYP_RUN_DIR"]).name,
                        r["endpoint"],
                        parent + ":" + r["strategy"] + ":" + r["subset"],
                        "-1",
                        r["model"] + ":" + metric,
                        r[metric],
                    ),
                )
                count += 1
write_json(out / "summary.json", dict(registered=count, config=cfg))
