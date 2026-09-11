"""Check interval-loss optimization against independently computed loss surfaces."""

import importlib.util
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location(
    "calibration_screen", Path(__file__).resolve().parents[1] / "scripts/calibration_screen.py"
)
cal = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cal)


def test_interval_shift_matches_dense_search_and_ignores_point_target():
    p = np.array([1.0, 2.0, 3.0, 4.0])
    lo, hi = p + np.array([0.2, 0.3, 0.5, 0.8]), p + np.array([0.4, 0.6, 0.7, 1.0])
    w = np.array([1.0, 2.0, 3.0, 1.0])
    model = cal.fit_map("shift_st", p, p + 10, lo, hi, w)
    value = np.dot(cal.soft(cal.apply_map(model, p), lo, hi), w)
    shifts = np.linspace(-1, 2, 3001)
    grid = np.maximum(lo[None, :] - (p[None, :] + shifts[:, None]), 0)
    grid += np.maximum((p[None, :] + shifts[:, None]) - hi[None, :], 0)
    assert value <= (grid @ w).min() + 1e-9
    assert model["a"] == 1.0


def test_interval_affine_finds_feasible_positive_map():
    p = np.arange(8, dtype=float)
    target = 0.7 * p + 0.9
    model = cal.fit_map("affine_st", p, target, target - 0.01, target + 0.01, np.ones(8))
    pred = cal.apply_map(model, p)
    assert np.all(np.diff(pred) > 0)
    assert np.max(np.abs(pred - target)) <= 0.01000001
