"""Gene-level differential expression for TNBC versus non-TNBC."""

from __future__ import annotations

import numpy as np
import pandas as pd

from tnbc_biomarkers.statistics import benjamini_hochberg, log2_fold_change, welch_test

TNBC = "TNBC"
NON_TNBC = "non-TNBC"


def differential_expression(
    expression: pd.DataFrame,
    groups: pd.Series,
    probe_for_gene: pd.Series | None = None,
    fdr_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
) -> pd.DataFrame:
    """Welch t-test of log2 expression, with Benjamini–Hochberg FDR.

    ``expression`` is samples by genes and must already be log2 intensity.
    The fold change is the TNBC mean minus the non-TNBC mean. A gene is called
    when FDR <= ``fdr_threshold`` and the absolute log2 fold change is at least
    ``log2fc_threshold``.
    """

    groups = groups.loc[expression.index]
    tnbc = expression.loc[groups == TNBC]
    other = expression.loc[groups == NON_TNBC]
    if tnbc.empty or other.empty:
        raise ValueError("Both TNBC and non-TNBC samples are required")
    log2fc = log2_fold_change(tnbc.to_numpy(), other.to_numpy())
    t_stat, p_value, df = welch_test(tnbc.to_numpy(), other.to_numpy())
    fdr = benjamini_hochberg(p_value)
    table = pd.DataFrame(
        {
            "gene": expression.columns.astype(str),
            "mean_tnbc": tnbc.mean(axis=0).to_numpy(),
            "mean_non_tnbc": other.mean(axis=0).to_numpy(),
            "log2_fold_change": log2fc,
            "t_statistic": t_stat,
            "welch_df": df,
            "p_value": p_value,
            "fdr": fdr,
        }
    )
    if probe_for_gene is not None:
        table.insert(1, "probe_id", probe_for_gene.reindex(table["gene"]).to_numpy())
    called = (
        np.isfinite(table["fdr"].to_numpy())
        & (table["fdr"].to_numpy() <= fdr_threshold)
        & (np.abs(table["log2_fold_change"].to_numpy()) >= log2fc_threshold)
    )
    table["significant"] = called
    table["higher_in"] = np.where(table["log2_fold_change"] >= 0, TNBC, NON_TNBC)
    table = table.sort_values(
        ["fdr", "p_value", "gene"],
        ascending=[True, True, True],
        na_position="last",
        kind="mergesort",
    ).reset_index(drop=True)
    return table


def passing_genes(table: pd.DataFrame, higher_in: str, limit: int) -> list[str]:
    """Genes called in one direction, ranked by FDR then absolute fold change."""

    subset = table.loc[table["significant"] & (table["higher_in"] == higher_in)].copy()
    subset["abs_log2_fold_change"] = subset["log2_fold_change"].abs()
    subset = subset.sort_values(
        ["fdr", "abs_log2_fold_change", "gene"],
        ascending=[True, False, True],
        kind="mergesort",
    )
    return subset["gene"].head(limit).tolist()
