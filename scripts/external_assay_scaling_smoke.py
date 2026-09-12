"""Execute real trainer task preparation and stop before model/GPU construction."""

import json
import os
from pathlib import Path

import train_external_assay as trainer


class PreparationComplete(Exception):
    pass


def stop_before_training(_seed):
    raise PreparationComplete


original_seed = trainer.torch.manual_seed
trainer.torch.manual_seed = stop_before_training
try:
    trainer.main()
except PreparationComplete:
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    tasks = json.loads((out / "tasks.json").read_text())
    assert tasks and all(t["scaling_policy"] == "pooled_shrinkage" for t in tasks)
    assert all(t["std"] >= 0.3 for t in tasks)
    assert (out / "training-variance-pools.parquet").exists()
    (out / "preparation-smoke.json").write_text(
        json.dumps(
            dict(
                passed=True,
                tasks=len(tasks),
                stage="Real trainer preparation complete; stopped before graph/model construction",
            ),
            indent=2,
        )
    )
finally:
    trainer.torch.manual_seed = original_seed
