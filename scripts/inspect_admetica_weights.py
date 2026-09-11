"""Inspect external CYP checkpoint metadata without executing pickled objects."""

import hashlib
import json
import os
import pickletools
import zipfile
from pathlib import Path

import torch

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
rows = []
for path in sorted(inp.glob("*.ckpt")):
    with zipfile.ZipFile(path) as z:
        data = z.read(next(n for n in z.namelist() if n.endswith("/data.pkl")))
    ops = list(pickletools.genops(data))
    strings = [
        arg for op, arg, pos in ops if op.name in ["BINUNICODE", "SHORT_BINUNICODE", "UNICODE"]
    ]
    primitive = []
    wanted = {
        "depth",
        "dropout",
        "d_h",
        "d_v",
        "d_e",
        "n_tasks",
        "input_dim",
        "hidden_dim",
        "batch_norm",
        "activation",
        "aggregation",
        "task_type",
        "dataset_type",
        "target_columns",
        "epochs",
        "epoch",
    }
    for i, (_op, arg, _pos) in enumerate(ops[:-1]):
        if isinstance(arg, str) and arg in wanted:
            following = [
                (o.name, a)
                for o, a, _ in ops[i + 1 : i + 6]
                if o.name not in ["BINPUT", "LONG_BINPUT", "MEMOIZE"]
            ]
            primitive.append({"key": arg, "following": following})
    row = dict(
        file=path.name,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        unsafe_globals=torch.serialization.get_unsafe_globals_in_checkpoint(path),
        metadata_fragments=primitive,
        tensor_keys=[s for s in strings if s.endswith((".weight", ".bias"))],
        relevant_strings=[
            s
            for s in strings
            if any(
                v in s.lower()
                for v in ["cyp", "challenge", "chemprop", "classification", "regression"]
            )
        ],
    )
    rows.append(row)
(out / "checkpoint-inspection.json").write_text(json.dumps(rows, indent=2, default=str))
print(json.dumps(rows, indent=2, default=str))
