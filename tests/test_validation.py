"""External transfer and replication rules."""

import numpy as np
import pandas as pd

from tnbc_biomarkers.config import AnalysisConfig
from tnbc_biomarkers.validation import replication_summary, replication_table, transfer_elastic_net
from tnbc_biomarkers.xena import ihc_group


def test_ihc_definition_excludes_equivocal_her2():
    assert ihc_group("Negative", "Negative", "Negative") == "TNBC"
    assert ihc_group("Positive", "Negative", "Negative") == "non-TNBC"
    assert ihc_group("Negative", "Negative", "Equivocal") is None
    assert ihc_group("Negative", "Indeterminate", "Negative") is None
    assert ihc_group("Negative", "Negative", None) is None


def test_replication_counts_direction_and_fdr():
    discovery = pd.DataFrame(
        {
            "gene": ["A", "B", "C", "D"],
            "log2_fold_change": [2.0, -1.5, 1.2, 1.1],
            "fdr": [0.01, 0.01, 0.01, 0.01],
            "higher_in": ["TNBC", "non-TNBC", "TNBC", "TNBC"],
            "significant": [True, True, True, True],
        }
    )
    external = pd.DataFrame(
        {
            "gene": ["A", "B", "C"],
            "log2_fold_change": [1.4, -0.4, -1.2],
            "fdr": [0.01, 0.2, 0.01],
            "higher_in": ["TNBC", "non-TNBC", "non-TNBC"],
            "significant": [True, False, True],
        }
    )
    table = replication_table(discovery, external, fdr_threshold=0.05, log2fc_threshold=1.0)
    summary = replication_summary(table)
    assert summary["n_discovery_called"] == 4
    assert summary["n_found_in_external"] == 3
    assert summary["n_same_direction"] == 2
    assert summary["n_same_direction_and_fdr"] == 1
    assert summary["n_same_direction_and_called"] == 1


def test_transfer_recovers_a_shared_signal():
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(6)]
    train = pd.DataFrame(rng.normal(scale=0.4, size=(40, 6)), columns=genes)
    test = pd.DataFrame(rng.normal(scale=0.4, size=(40, 6)), columns=genes)
    train.iloc[:20, 0] += 3
    train.iloc[20:, 0] -= 3
    test.iloc[:20, 0] += 3
    test.iloc[20:, 0] -= 3
    groups = pd.Series(["TNBC"] * 20 + ["non-TNBC"] * 20)
    config = AnalysisConfig(
        outer_splits=2,
        inner_splits=2,
        n_features=3,
        cs=(0.1, 1.0),
        l1_ratios=(0.5,),
        random_state=0,
        min_samples_per_group=4,
    )
    result = transfer_elastic_net(train, groups, test, groups, config)
    assert result["elastic_net"]["roc_auc"] > 0.8
    assert result["external_majority_baseline"]["balanced_accuracy"] == 0.5
    assert "G0" in result["genes"]
