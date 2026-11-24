"""Sample correlation filtering."""

import numpy as np
import pandas as pd

from tnbc_biomarkers.qc import sample_quality, value_scale_summary


def test_low_correlation_sample_is_flagged():
    rng = np.random.default_rng(2)
    shared = rng.normal(size=(1, 30))
    values = np.repeat(shared, 6, axis=0) + rng.normal(scale=0.05, size=(6, 30))
    values[5] = rng.normal(size=30)
    expression = pd.DataFrame(values, index=[f"S{i}" for i in range(6)])
    quality = sample_quality(expression, correlation_threshold=0.80)
    assert quality["pass_correlation"].tolist() == [True, True, True, True, True, False]
    assert value_scale_summary(expression)["assumed_log2"] is True


def test_linear_scale_is_not_treated_as_log2():
    expression = pd.DataFrame(np.full((4, 3), 100.0))
    summary = value_scale_summary(expression)
    assert summary["assumed_log2"] is False
    assert summary["median"] == 100.0
