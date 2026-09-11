"""Connectivity-aware exclusion for external-data regression comparisons."""

import numpy as np


def snapshot_input(directory, filename):
    """Match the original basename exactly after the run's numeric prefix."""
    matches = [p for p in directory.iterdir() if p.name.partition("-")[2] == filename]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one input named {filename}; found {len(matches)}")
    return matches[0]


def grouped_membership(
    ids,
    connectivity,
    folds,
    outer_fold,
    seed,
    fraction,
    outer_connectivity=(),
    test_connectivity=(),
):
    """Preserve fixed outer folds; exclude their connectivity from all fitting roles."""
    ids, keys, folds = np.asarray(ids), np.asarray(connectivity), np.asarray(folds)
    if len(set(ids)) != len(ids) or not 0 < fraction < 1:
        raise ValueError("Unique molecular IDs and a validation fraction in (0,1) required")
    held = np.flatnonzero(folds == outer_fold)
    forbidden = set(keys[held]) | set(outer_connectivity) | set(test_connectivity)
    eligible_keys = sorted(set(keys) - forbidden)
    shuffled = np.random.default_rng(seed).permutation(eligible_keys)
    n_val = max(1, int(len(shuffled) * fraction))
    val_keys = set(shuffled[:n_val])
    train_keys = set(shuffled[n_val:])
    train = np.flatnonzero(np.isin(keys, list(train_keys)))
    val = np.flatnonzero(np.isin(keys, list(val_keys)))
    if not len(train) or not len(val) or not len(held):
        raise ValueError("Empty train, validation or held-out partition")
    assert set(keys[train]).isdisjoint(keys[val])
    assert (set(keys[train]) | set(keys[val])).isdisjoint(forbidden)
    return train, val, held
