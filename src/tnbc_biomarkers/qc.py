"""Sample-level quality control for a gene-expression matrix."""

from __future__ import annotations

import numpy as np
import pandas as pd


def sample_quality(expression: pd.DataFrame, correlation_threshold: float) -> pd.DataFrame:
    """Summarise each sample and flag low median pairwise correlation.

    ``expression`` is samples by genes. A sample fails when its median Pearson
    correlation with the other samples is below ``correlation_threshold`` or is
    undefined.
    """

    values = expression.to_numpy(dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        correlation = np.corrcoef(values)
    np.fill_diagonal(correlation, np.nan)
    median_correlation = np.nanmedian(correlation, axis=1)
    frame = pd.DataFrame(
        {
            "median_expression": np.nanmedian(values, axis=1),
            "iqr": np.nanpercentile(values, 75, axis=1) - np.nanpercentile(values, 25, axis=1),
            "median_correlation": median_correlation,
        },
        index=expression.index,
    )
    frame["pass_correlation"] = frame["median_correlation"] >= correlation_threshold
    return frame


def value_scale_summary(expression: pd.DataFrame) -> dict[str, float | str]:
    """Describe the matrix scale. Submitted RMA values are left unchanged.

    A median above 30 would indicate a linear scale and is reported so the
    caller can refuse to treat those values as log2 intensities.
    """

    values = expression.to_numpy(dtype=float)
    median = float(np.nanmedian(values))
    return {
        "min": float(np.nanmin(values)),
        "median": median,
        "max": float(np.nanmax(values)),
        "assumed_log2": bool(median <= 30),
    }
