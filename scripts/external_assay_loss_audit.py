"""Independently check ranking orientation and interval-loss gradients on censored cases."""

import json
import os
from pathlib import Path

import numpy as np
import torch
from external_loss_utils import interval_mse, ordered_pairs, ranking_loss

rng = np.random.default_rng(83)
checked = 0
for n in [0, 1, 2, 17, 1000]:
    lower = rng.uniform(2, 9, n)
    upper = lower + rng.uniform(0, 2, n)
    lower[rng.random(n) < 0.2] = -np.inf
    upper[rng.random(n) < 0.2] = np.inf
    pairs = ordered_pairs(lower, upper, 64, 0.25, rng)
    assert len(pairs) <= 64
    for strong, weak in pairs:
        assert strong != weak and lower[strong] >= upper[weak] + 0.25
    checked += len(pairs)
p = torch.tensor([4.0, 8.0, 6.0, 5.0, 7.0], requires_grad=True)
lo = torch.tensor([5.0, -float("inf"), 5.0, 5.0, -float("inf")])
hi = torch.tensor([float("inf"), 7.0, 7.0, 5.0, float("inf")])
loss = interval_mse(p, lo, hi)
loss.backward()
assert torch.isfinite(loss) and torch.isfinite(p.grad).all()
assert torch.allclose(p.grad, torch.tensor([-0.4, 0.4, 0.0, 0.0, 0.0]))
s = torch.tensor([7.0], requires_grad=True)
w = torch.tensor([5.0], requires_grad=True)
ranking_loss(s, w).backward()
assert s.grad.item() < 0 and w.grad.item() > 0
assert ranking_loss(torch.tensor([8.0]), torch.tensor([5.0])) < ranking_loss(
    torch.tensor([7.0]), torch.tensor([5.0])
)
result = dict(
    pair_constraints_checked=checked,
    interval_loss=float(loss.detach()),
    interval_gradients=p.grad.tolist(),
    ranking_gradient_direction="increases stronger prediction and decreases weaker",
    passed=True,
)
Path(os.environ["CYP_OUTPUT_DIR"], "loss-audit.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result))
