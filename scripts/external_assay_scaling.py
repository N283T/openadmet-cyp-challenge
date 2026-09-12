"""Training-only variance priors for optional assay scale shrinkage."""

import numpy as np
import polars as pl

GROUP_COLUMNS = ["source", "isoform", "measurement_type"]


def training_variance_pools(observations, train_compound_ids, minimum_std=0.3):
    """Pool within-assay residuals after excluding every validation compound.

    Input is already deduplicated per compound/assay using the trainer policy.
    Assay offsets therefore do not inflate the residual variance prior.
    """
    if minimum_std <= 0:
        raise ValueError("minimum_std must be positive")
    exact = observations.filter(
        pl.col("compound_id").is_in(pl.Series(train_compound_ids).implode())
        & pl.col("lower_bound").is_finite()
        & (pl.col("lower_bound") == pl.col("upper_bound"))
    )
    exact = exact.with_columns(
        (pl.col("lower_bound") - pl.col("lower_bound").mean().over("assay_id")).alias(
            "assay_residual"
        )
    )
    return (
        exact.group_by(GROUP_COLUMNS)
        .agg(
            pl.len().alias("pool_exact_n"),
            pl.col("assay_id").n_unique().alias("pool_assays"),
            (pl.col("assay_residual") ** 2).mean().alias("pool_raw_variance"),
        )
        .with_columns(
            pl.col("pool_raw_variance").clip(lower_bound=minimum_std**2).alias("pool_variance")
        )
        .sort(GROUP_COLUMNS)
    )


def shrunk_assay_std(exact_training_values, pool_variance, pseudocount=20.0, minimum_std=0.3):
    """Shrink population variance with a fixed training-only prior, retaining assay mean."""
    if pseudocount <= 0 or minimum_std <= 0:
        raise ValueError("pseudocount and minimum_std must be positive")
    if not np.isfinite(pool_variance) or pool_variance < 0:
        raise ValueError("pool_variance must be finite and nonnegative")
    values = np.asarray(exact_training_values, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Exact training values must be finite")
    n = len(values)
    variance = float(np.var(values)) if n else 0.0
    result = (n * variance + pseudocount * pool_variance) / (n + pseudocount)
    return max(minimum_std, float(np.sqrt(result)))
