"""Record whether a pinned GLORYxR checkout is usable for model evaluation."""

import json
import os
from pathlib import Path

import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
tutorial = json.loads((inp / cfg["tutorial"]).read_text())
tutorial_text = "\n".join(
    line for cell in tutorial["cells"] for line in cell.get("source", []) if isinstance(line, str)
)
rule_rows = []
for filename in cfg["reaction_rules"]:
    table = pl.read_csv(inp / filename, infer_schema_length=10000, truncate_ragged_lines=True)
    rule_rows.append({"file": filename, "rows": table.height, "columns": len(table.columns)})
pl.DataFrame(rule_rows).write_csv(out / "reaction-rule-inventory.csv")
full_models_published = "Zenodo archive (to be published)" not in tutorial_text
example_warning = "NOT trained on the full GLORYxR dataset" in tutorial_text
summary = {
    "source_commit": cfg["source_commit"],
    "python_requirement": ">=3.13",
    "reaction_rule_files": rule_rows,
    "bundled_joblib_models": 1,
    "bundled_model_is_tutorial_only": example_warning,
    "full_model_archive_available_in_pinned_source": full_models_published,
    "decision": "blocked_model_evaluation"
    if example_warning and not full_models_published
    else "ready",
    "reason": (
        "The only bundled model is explicitly a rudimentary tutorial model; "
        "the full-model Zenodo archive is still marked to be published."
    ),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
