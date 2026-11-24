"""Checks for fold change, the Welch test, and Benjamini–Hochberg FDR."""

import numpy as np
import pandas as pd
from scipy import stats

from tnbc_biomarkers.differential import differential_expression
from tnbc_biomarkers.statistics import benjamini_hochberg, log2_fold_change, welch_test


def test_benjamini_hochberg_known_values():
    p_values = np.array([0.01, 0.04, 0.03])
    adjusted = benjamini_hochberg(p_values)
    assert adjusted == pytest_approx([0.03, 0.04, 0.04])


def test_benjamini_hochberg_ignores_non_finite_and_is_aligned():
    p_values = np.array([np.nan, 0.2, 0.01])
    adjusted = benjamini_hochberg(p_values)
    assert np.isnan(adjusted[0])
    assert adjusted[1] == pytest_approx(0.2)
    assert adjusted[2] == pytest_approx(0.02)
    assert adjusted[2] <= adjusted[1]


def test_welch_matches_hand_formula():
    group_a = np.array([[1.0], [2.0], [3.0]])
    group_b = np.array([[4.0], [6.0], [8.0]])
    t_stat, p_value, df = welch_test(group_a, group_b)
    variance_a = 1.0
    variance_b = 4.0
    se = np.sqrt(variance_a / 3 + variance_b / 3)
    expected_t = (2.0 - 6.0) / se
    expected_df = (variance_a / 3 + variance_b / 3) ** 2 / (
        (variance_a / 3) ** 2 / 2 + (variance_b / 3) ** 2 / 2
    )
    expected_p = 2 * stats.t.sf(abs(expected_t), expected_df)
    assert t_stat[0] == pytest_approx(expected_t)
    assert df[0] == pytest_approx(expected_df)
    assert p_value[0] == pytest_approx(expected_p)
    assert log2_fold_change(group_a, group_b)[0] == pytest_approx(-4.0)


def test_contrast_is_tnbc_minus_non_tnbc():
    expression = pd.DataFrame(
        {"G1": [1.0, 1.0, 5.0, 5.0], "G2": [2.0, 3.0, 2.5, 2.5]},
        index=["s1", "s2", "s3", "s4"],
    )
    groups = pd.Series(["TNBC", "TNBC", "non-TNBC", "non-TNBC"], index=expression.index)
    table = differential_expression(expression, groups).set_index("gene")
    assert table.loc["G1", "log2_fold_change"] == pytest_approx(-4.0)
    assert table.loc["G1", "higher_in"] == "non-TNBC"
    assert table.loc["G1", "significant"]
    assert table.loc["G2", "log2_fold_change"] == pytest_approx(0.0)
    assert not table.loc["G2", "significant"]


def pytest_approx(values):
    import pytest

    return pytest.approx(values)
