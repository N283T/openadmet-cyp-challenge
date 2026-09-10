import pytest
import torch

from cyp_challenge.neural_losses import masked_loss


def test_missing_labels_have_zero_gradient_and_heads_have_equal_weight():
    pred = torch.tensor([[1.0, 2.0, 100.0], [3.0, 999.0, 200.0]], requires_grad=True)
    target = torch.tensor([[0.0, 0.0, float("nan")], [0.0, float("nan"), float("nan")]])
    loss = masked_loss(pred, target, False)
    assert float(loss.detach()) == pytest.approx(((1 + 9) / 2 + 4) / 2)
    loss.backward()
    assert torch.isfinite(pred.grad).all()
    assert pred.grad[1, 1] == 0
    assert (pred.grad[:, 2] == 0).all()


def test_binary_loss_masks_missing_heads_and_rejects_empty_batch():
    pred = torch.zeros((2, 2), requires_grad=True)
    target = torch.tensor([[1.0, float("nan")], [0.0, float("nan")]])
    loss = masked_loss(pred, target, True)
    assert float(loss.detach()) == pytest.approx(0.69314718)
    loss.backward()
    assert (pred.grad[:, 1] == 0).all()
    with pytest.raises(ValueError):
        masked_loss(pred, torch.full((2, 2), float("nan")), True)
