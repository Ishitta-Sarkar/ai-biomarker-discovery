"""Fixed analysis choices for the GSE76275 TNBC contrast.

These defaults are the methods reported in the README. Tests pass a modified
``AnalysisConfig`` so the continuous-integration run stays small.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class AnalysisConfig:
    """Parameters for quality control, differential expression, and models."""

    accession: str = "GSE76275"
    fdr_threshold: float = 0.05
    log2fc_threshold: float = 1.0
    correlation_threshold: float = 0.80
    pca_genes: int = 2000
    heatmap_genes: int = 40
    volcano_labels: int = 10
    outer_splits: int = 5
    inner_splits: int = 5
    n_features: int = 250
    enrichment_cap: int = 500
    random_state: int = 42
    cs: tuple[float, ...] = (1e-2, 1e-1, 1.0, 10.0, 100.0)
    l1_ratios: tuple[float, ...] = (0.15, 0.5, 0.85)
    enrichment_sources: tuple[str, ...] = ("GO:BP", "KEGG", "REAC")
    nonzero_coef: float = 1e-8
    min_samples_per_group: int = 20
