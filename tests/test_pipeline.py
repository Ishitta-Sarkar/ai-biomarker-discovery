"""End-to-end run on a synthetic matrix, with enrichment replaced by a stub."""

from pathlib import Path

import numpy as np
import pandas as pd

from tnbc_biomarkers.config import AnalysisConfig
from tnbc_biomarkers.pipeline import analyse_expression


def _enrich(genes, background, sources):
    table = pd.DataFrame(
        [
            {
                "source": "GO:BP",
                "native": "GO:0008150",
                "name": "biological process",
                "p_value": 0.01,
                "term_size": 10,
                "query_size": len(genes),
                "intersection_size": min(2, len(genes)),
                "precision": 0.5,
                "recall": 0.1,
            }
        ]
    )
    meta = {
        "version": "stub",
        "timestamp": "2026-11-24T00:00:00+00:00",
        "domain_scope": "custom",
        "background_fallback": False,
    }
    if not genes:
        table = table.iloc[0:0]
    return table, meta


def test_analyse_expression_writes_figures_and_tables(tmp_path: Path):
    rng = np.random.default_rng(0)
    probes = [f"p{i}" for i in range(20)]
    sample_ids = [f"S{i:02d}" for i in range(24)]
    values = rng.normal(loc=6.0, scale=0.3, size=(20, 24))
    values += rng.normal(size=(20, 1))
    values[0, :12] += 2.0
    values[0, 12:] -= 2.0
    expression = pd.DataFrame(values, index=probes, columns=sample_ids)
    symbols = pd.Series({probe: f"G{i}" for i, probe in enumerate(probes)})
    samples = pd.DataFrame(
        {
            "title": sample_ids,
            "series": ["GSE1"] * 12 + ["GSE2"] * 12,
            "group": ["TNBC"] * 12 + ["non-TNBC"] * 12,
            "characteristics": [{}] * 24,
            "tn_status_raw": ["TN"] * 12 + ["not TN"] * 12,
            "er": ["Negative"] * 12 + ["Positive"] * 12,
            "pr": ["Negative"] * 24,
            "her2": ["Negative"] * 24,
            "receptor_group": ["TNBC"] * 12 + ["non-TNBC"] * 12,
            "data_processing": ["RMA"] * 24,
        },
        index=sample_ids,
    )
    config = AnalysisConfig(
        outer_splits=2,
        inner_splits=2,
        n_features=5,
        pca_genes=10,
        heatmap_genes=6,
        enrichment_cap=5,
        correlation_threshold=0.0,
        cs=(0.1, 1.0),
        l1_ratios=(0.5,),
        random_state=0,
        min_samples_per_group=6,
    )
    summary = analyse_expression(
        expression,
        samples,
        symbols,
        tmp_path,
        config,
        enrich_fn=_enrich,
    )
    assert (tmp_path / "figures" / "volcano.png").stat().st_size > 0
    assert (tmp_path / "figures" / "heatmap_top_genes.png").stat().st_size > 0
    assert (tmp_path / "figures" / "pca_samples.png").stat().st_size > 0
    assert summary["qc"]["n_tnbc"] == 12
    assert summary["differential_expression"]["n_genes_tested"] == 20
    assert summary["classification"]["elastic_net"]["pooled_roc_auc"] > 0.8
    saved = pd.read_csv(tmp_path / "tables" / "differential_expression.csv")
    assert {"gene", "log2_fold_change", "fdr", "significant"}.issubset(saved.columns)
