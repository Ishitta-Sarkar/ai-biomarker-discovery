"""Nested elastic-net classification on a separable synthetic cohort."""

import numpy as np
import pandas as pd

from tnbc_biomarkers.classify import nested_classification, select_feature_indices
from tnbc_biomarkers.config import AnalysisConfig


def test_feature_selection_uses_only_the_training_rows():
    x_train = np.array(
        [
            [0.0, 5.0],
            [0.1, 5.2],
            [3.0, 0.0],
            [3.2, 0.1],
        ]
    )
    y_train = np.array([1, 1, 0, 0])
    selected = select_feature_indices(x_train, y_train, k=1)
    assert selected.tolist() == [1]


def test_nested_elastic_net_beats_the_majority_baseline():
    rng = np.random.default_rng(1)
    n_samples = 40
    x = rng.normal(scale=0.4, size=(n_samples, 12))
    x[:20, 0] += 3
    x[20:, 0] -= 3
    expression = pd.DataFrame(x, columns=[f"G{i}" for i in range(12)])
    groups = pd.Series(["TNBC"] * 20 + ["non-TNBC"] * 20)
    config = AnalysisConfig(
        outer_splits=2,
        inner_splits=2,
        n_features=4,
        cs=(0.1, 1.0),
        l1_ratios=(0.5,),
        random_state=1,
        min_samples_per_group=2,
    )
    result = nested_classification(expression, groups, config)
    assert result["elastic_net"]["pooled_roc_auc"] > 0.8
    assert result["elastic_net"]["pooled_roc_auc"] > result["majority_baseline"]["pooled_roc_auc"]
    assert result["majority_baseline"]["mean"]["balanced_accuracy"] == pytest_approx(0.5)
    assert 0.0 <= result["l2_logistic"]["pooled_roc_auc"] <= 1.0
    assert result["convergence_warnings"] == 0


def pytest_approx(value):
    import pytest

    return pytest.approx(value)
