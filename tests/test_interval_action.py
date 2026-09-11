import importlib.util
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location(
    "interval_readout", Path(__file__).resolve().parents[1] / "scripts/interval_readout.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_marginal_mixture_median_minimizes_expected_interval_distance():
    lo = np.array([[0.0, 1.0, 2.0, 7.0], [-2.0, -1.0, 0.0, 1.0]])
    hi = np.array([[0.5, 3.0, 4.0, 8.0], [3.0, 4.0, 5.0, 6.0]])
    action = module.mixture_median(lo, hi)
    grid = np.linspace(-3, 9, 12001)
    for i in range(len(lo)):
        errors = (
            np.maximum(lo[i][None, :] - grid[:, None], 0)
            + np.maximum(grid[:, None] - hi[i][None, :], 0)
        ).mean(axis=1)
        actual = (np.maximum(lo[i] - action[i], 0) + np.maximum(action[i] - hi[i], 0)).mean()
        assert actual <= errors.min() + 1e-10


def test_half_weight_endpoint_mae_equals_interval_loss_plus_constant():
    lo = np.array([-2.0, 1.0, 3.0])
    hi = np.array([0.0, 2.0, 7.0])
    for pred in np.linspace(-3, 8, 25):
        endpoint_mae = 0.5 * (np.abs(pred - lo) + np.abs(pred - hi)).mean()
        interval_loss = (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).mean()
        assert np.isclose(endpoint_mae - interval_loss, 0.5 * (hi - lo).mean())
