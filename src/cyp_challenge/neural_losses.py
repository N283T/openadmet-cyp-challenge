"""Masked multitask objectives with equal contribution from observed heads."""

import torch


def masked_loss(prediction, target, binary):
    mask = torch.isfinite(target)
    truth = torch.nan_to_num(target)
    loss = (
        torch.nn.functional.binary_cross_entropy_with_logits(prediction, truth, reduction="none")
        if binary
        else (prediction - truth) ** 2
    )
    counts = mask.sum(dim=0)
    present = counts > 0
    if not bool(present.any()):
        raise ValueError("No observed labels in batch")
    return ((loss * mask).sum(dim=0)[present] / counts[present]).mean()
