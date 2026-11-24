"""Two-group statistics used by differential expression and the classifier."""

from __future__ import annotations

import numpy as np
from scipy import stats


def log2_fold_change(group_a: np.ndarray, group_b: np.ndarray) -> np.ndarray:
    """Difference of column means.

    Both matrices must already be on a log2 scale. The result is the log2 fold
    change of group A relative to group B (positive means higher in A).
    """

    return np.asarray(group_a, dtype=float).mean(axis=0) - np.asarray(group_b, dtype=float).mean(
        axis=0
    )


def welch_df(group_a: np.ndarray, group_b: np.ndarray) -> np.ndarray:
    """Welch–Satterthwaite degrees of freedom for each column."""

    a = np.asarray(group_a, dtype=float)
    b = np.asarray(group_b, dtype=float)
    n1 = a.shape[0]
    n2 = b.shape[0]
    v1 = a.var(axis=0, ddof=1)
    v2 = b.var(axis=0, ddof=1)
    se2 = v1 / n1 + v2 / n2
    numerator = se2**2
    denominator = (v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        df = numerator / denominator
    return df


def welch_test(
    group_a: np.ndarray, group_b: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Two-sided Welch t-test computed independently for each column.

    Returns the t statistic, the two-sided p-value, and the Welch–Satterthwaite
    degrees of freedom. Columns with undefined variance yield non-finite values.
    """

    a = np.asarray(group_a, dtype=float)
    b = np.asarray(group_b, dtype=float)
    t_stat, p_value = stats.ttest_ind(a, b, axis=0, equal_var=False, alternative="two-sided")
    return np.asarray(t_stat, dtype=float), np.asarray(p_value, dtype=float), welch_df(a, b)


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR adjusted p-values, preserving the input order.

    Non-finite p-values stay missing and are excluded from the correction.
    Adjusted values are capped at 1.
    """

    p = np.asarray(p_values, dtype=float)
    adjusted = np.full(p.shape, np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return adjusted
    pv = p[finite]
    n = pv.size
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order]
    running = 1.0
    adj_ranked = np.empty(n, dtype=float)
    for i in range(n - 1, -1, -1):
        running = min(running, ranked[i] * n / (i + 1))
        adj_ranked[i] = running
    restored = np.empty(n, dtype=float)
    restored[order] = np.clip(adj_ranked, 0.0, 1.0)
    adjusted[finite] = restored
    return adjusted
