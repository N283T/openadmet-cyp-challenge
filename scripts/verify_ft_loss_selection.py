"""Validate FT selection scores, saved epochs, default loss and holdout separation."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from finetune_external_multitask import inner_pic50_score, masked_loss

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
p = torch.tensor([[1.0, 2.0], [3.0, 4.0]])
y = torch.tensor([[0.0, float("nan")], [2.0, 6.0]])
w = torch.tensor([1.0, 2.0])
mask = torch.isfinite(y)
old = ((p - torch.nan_to_num(y)).square() * mask).sum(0) / mask.sum(0).clamp_min(1)
old = (old[mask.any(0)] * w[mask.any(0)]).sum() / w[mask.any(0)].sum()
assert torch.equal(old, masked_loss(p, y, w))
assert torch.isclose(masked_loss(p, y, w, "mae"), torch.tensor(5 / 3))
truth = np.array([[1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 4.0, 5.0], [3.0, 4.0, 5.0, 6.0]])
pred = truth + np.array([0.1, 0.2, 0.3, 0.4])
lo, hi = truth - 0.15, truth + 0.15
score, perhead = inner_pic50_score(scorer, pred, truth, lo, hi)
manual = [
    scorer.rae_soft_threshold_absolute_error(
        truth[:, i], pred[:, i], y_true_lower=lo[:, i], y_true_upper=hi[:, i]
    )
    for i in range(4)
]
assert score == float(np.mean(manual)) and perhead == manual
checks = []
for item in cfg["folds"]:
    history = [json.loads(line) for line in (inp / item["history"]).read_text().splitlines()]
    assert len(history) == 2
    for row in history:
        assert np.isfinite(row["validation_mse"]) and np.isfinite(row["validation_pic50_st_rae"])
        assert np.isclose(
            row["validation_pic50_st_rae"],
            np.mean(list(row["validation_pic50_head_st_rae"].values())),
        )
        assert row["selection_score"] == row["validation_pic50_st_rae"]
    best = min(history, key=lambda r: r["selection_score"])
    checkpoint = torch.load(inp / item["checkpoint"], map_location="cpu", weights_only=False)
    assert checkpoint["epoch"] == best["epoch"]
    members = pl.read_parquet(inp / item["membership"])
    roles = set(members["role"])
    print("roles", roles, flush=True)
    held = members.filter(pl.col("role") == "outer_holdout")
    iv = members.filter(pl.col("role") == "inner_validation")
    assert iv.height > 0 and held.height > 0
    assert set(held["compound_id"]).isdisjoint(set(iv["compound_id"]))
    checks.append(
        dict(
            fold=item["fold"],
            selected_epoch=best["epoch"],
            selection=best["selection_score"],
            inner_rows=iv.height,
            outer_rows=held.height,
        )
    )
(out / "summary.json").write_text(
    json.dumps(
        dict(
            default_mse_bitwise_equal=True,
            mae_formula=True,
            official_macro_formula=True,
            fold_checks=checks,
        ),
        indent=2,
    )
)
