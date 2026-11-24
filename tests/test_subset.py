"""The cached GSE76275 subset exercises the same statistics as the full run."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tnbc_biomarkers.classify import nested_classification
from tnbc_biomarkers.config import AnalysisConfig
from tnbc_biomarkers.differential import differential_expression

FIXTURE = Path(__file__).parent / "fixtures" / "gse76275_subset.csv"
PROVENANCE = FIXTURE.with_suffix(".provenance.json")


def test_cached_subset_matches_direct_fold_changes():
    table = pd.read_csv(FIXTURE)
    provenance = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    assert provenance["source_accession"] == "GSE76275"
    assert provenance["platform"] == "GPL570"
    assert set(table["group"]) == {"TNBC", "non-TNBC"}
    expression = table.drop(columns=["sample_id", "group"])
    values = expression.to_numpy(dtype=float)
    assert np.any(np.abs(values - np.round(values)) > 1e-6)
    groups = pd.Series(table["group"].to_numpy(), index=table["sample_id"])
    expression.index = table["sample_id"]
    result = differential_expression(expression, groups).set_index("gene")
    for gene in expression.columns:
        tnbc = expression.loc[groups.to_numpy() == "TNBC", gene].mean()
        other = expression.loc[groups.to_numpy() == "non-TNBC", gene].mean()
        assert result.loc[gene, "log2_fold_change"] == pytest.approx(tnbc - other)
    assert result["fdr"].notna().any()
    config = AnalysisConfig(
        outer_splits=2,
        inner_splits=2,
        n_features=5,
        cs=(0.1, 1.0),
        l1_ratios=(0.5,),
        random_state=0,
        min_samples_per_group=2,
    )
    fitted = nested_classification(expression, groups, config)
    for name in ("elastic_net", "l2_logistic", "majority_baseline"):
        auc = fitted[name]["pooled_roc_auc"]
        assert 0.0 <= auc <= 1.0
        assert 0.0 <= fitted[name]["mean"]["balanced_accuracy"] <= 1.0
