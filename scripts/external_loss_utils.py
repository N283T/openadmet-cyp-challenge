"""Assay-local interval and provably ordered pair losses."""

import numpy as np
import torch
import torch.nn.functional as F


def interval_mse(prediction, lower, upper):
    return (torch.relu(lower - prediction) ** 2 + torch.relu(prediction - upper) ** 2).mean()


def ordered_pairs(lower, upper, maximum, gap, rng):
    n = len(lower)
    if n < 2:
        return np.empty((0, 2), dtype=np.int64)
    if n * (n - 1) // 2 <= maximum * 8:
        a, b = np.triu_indices(n, 1)
    else:
        a, b = rng.integers(0, n, size=(2, maximum * 20))
    forward = (a != b) & (lower[a] - upper[b] >= gap)
    reverse = (a != b) & (lower[b] - upper[a] >= gap)
    pairs = np.unique(
        np.concatenate(
            [np.column_stack([a[forward], b[forward]]), np.column_stack([b[reverse], a[reverse]])]
        ),
        axis=0,
    )
    if len(pairs) > maximum:
        pairs = pairs[rng.choice(len(pairs), maximum, replace=False)]
    return pairs


def ranking_loss(stronger, weaker):
    return F.softplus(-(stronger - weaker)).mean()
