"""External test of a GSE76275 model and replication of its called genes."""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from tnbc_biomarkers.classify import (
    _logistic,
    _scalar,
    classification_metrics,
    nested_classification,
    select_feature_indices,
)
from tnbc_biomarkers.config import AnalysisConfig
from tnbc_biomarkers.differential import differential_expression
from tnbc_biomarkers.plots import plot_roc
from tnbc_biomarkers.qc import sample_quality, value_scale_summary

TNBC = "TNBC"


def transfer_elastic_net(
    train: pd.DataFrame,
    train_groups: pd.Series,
    test: pd.DataFrame,
    test_groups: pd.Series,
    config: AnalysisConfig,
) -> dict[str, object]:
    """Fit an elastic net on the discovery cohort and score the external cohort.

    Genes are ranked by the absolute Welch t statistic on the discovery cohort
    only. Each cohort is z-scored on its own samples before the coefficients,
    which were learned in discovery z-score units, are applied to the external
    cohort. The external labels are not used to choose genes or hyperparameters.
    """

    shared = train.columns.intersection(test.columns)
    if len(shared) < 2:
        raise ValueError("The cohorts share fewer than two genes")
    y_train = (train_groups.loc[train.index] == TNBC).astype(int).to_numpy()
    y_test = (test_groups.loc[test.index] == TNBC).astype(int).to_numpy()
    training_matrix = train.loc[:, shared].to_numpy(dtype=float)
    selected = select_feature_indices(training_matrix, y_train, config.n_features)
    genes = shared.to_numpy()[selected]
    x_train = StandardScaler().fit_transform(train.loc[:, genes].to_numpy(dtype=float))
    x_test = StandardScaler().fit_transform(test.loc[:, genes].to_numpy(dtype=float))
    inner = StratifiedKFold(
        n_splits=config.inner_splits,
        shuffle=True,
        random_state=config.random_state,
    )
    model = _logistic(config, elastic_net=True, inner=inner)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        model.fit(x_train, y_train)
        n_warnings = sum(issubclass(item.category, ConvergenceWarning) for item in caught)
    scores = model.predict_proba(x_test)[:, 1]
    coefficients = np.asarray(model.coef_, dtype=float).reshape(-1)
    metrics = classification_metrics(y_test, scores)
    training_majority = int(np.bincount(y_train).argmax())
    external_majority = int(np.bincount(y_test).argmax())
    training_baseline = classification_metrics(
        y_test,
        np.full(y_test.shape, float(y_train.mean())),
        np.full(y_test.shape, training_majority),
    )
    external_baseline = classification_metrics(
        y_test,
        np.full(y_test.shape, float(y_test.mean())),
        np.full(y_test.shape, external_majority),
    )
    return {
        "n_shared_genes": int(len(shared)),
        "n_features": int(len(genes)),
        "genes": [str(gene) for gene in genes],
        "C": _scalar(model.C_),
        "l1_ratio": _scalar(model.l1_ratio_),
        "n_nonzero": int(np.sum(np.abs(coefficients) > config.nonzero_coef)),
        "convergence_warnings": n_warnings,
        "scaling": "Each gene was z-scored inside its own cohort.",
        "positive_class": TNBC,
        "decision_threshold": 0.5,
        "n_test_samples": int(len(y_test)),
        "n_test_tnbc": int(y_test.sum()),
        "n_test_non_tnbc": int(len(y_test) - y_test.sum()),
        "elastic_net": metrics,
        "training_majority_baseline": training_baseline,
        "external_majority_baseline": external_baseline,
        "y_true": y_test,
        "scores": scores,
    }


def replication_table(
    discovery: pd.DataFrame,
    external: pd.DataFrame,
    fdr_threshold: float,
    log2fc_threshold: float,
) -> pd.DataFrame:
    """Compare each discovery call with the external log2 fold change and FDR."""

    called = discovery.loc[discovery["significant"]].copy()
    external_by_gene = external.drop_duplicates("gene", keep="first").set_index("gene")
    rows = []
    for record in called.itertuples(index=False):
        gene = record.gene
        if gene not in external_by_gene.index:
            rows.append(
                {
                    "gene": gene,
                    "discovery_log2_fold_change": record.log2_fold_change,
                    "discovery_fdr": record.fdr,
                    "discovery_higher_in": record.higher_in,
                    "found_in_external": False,
                    "external_log2_fold_change": np.nan,
                    "external_fdr": np.nan,
                    "same_direction": False,
                    "external_fdr_significant": False,
                    "external_called": False,
                }
            )
            continue
        other = external_by_gene.loc[gene]
        same_direction = bool(
            np.sign(record.log2_fold_change) == np.sign(other["log2_fold_change"])
            and np.sign(record.log2_fold_change) != 0
        )
        fdr_significant = bool(other["fdr"] <= fdr_threshold)
        external_called = bool(
            fdr_significant and abs(other["log2_fold_change"]) >= log2fc_threshold
        )
        rows.append(
            {
                "gene": gene,
                "discovery_log2_fold_change": record.log2_fold_change,
                "discovery_fdr": record.fdr,
                "discovery_higher_in": record.higher_in,
                "found_in_external": True,
                "external_log2_fold_change": float(other["log2_fold_change"]),
                "external_fdr": float(other["fdr"]),
                "same_direction": same_direction,
                "external_fdr_significant": fdr_significant,
                "external_called": external_called,
            }
        )
    table = pd.DataFrame(rows)
    return table


def replication_summary(table: pd.DataFrame) -> dict[str, int]:
    if table.empty:
        return {
            "n_discovery_called": 0,
            "n_found_in_external": 0,
            "n_same_direction": 0,
            "n_same_direction_and_fdr": 0,
            "n_same_direction_and_called": 0,
        }
    found = table["found_in_external"]
    return {
        "n_discovery_called": int(table.shape[0]),
        "n_found_in_external": int(found.sum()),
        "n_same_direction": int((found & table["same_direction"]).sum()),
        "n_same_direction_and_fdr": int(
            (found & table["same_direction"] & table["external_fdr_significant"]).sum()
        ),
        "n_same_direction_and_called": int(
            (found & table["same_direction"] & table["external_called"]).sum()
        ),
    }


def within_batch_summary(result: dict) -> dict:
    """JSON-ready nested-CV block without coefficient frames."""

    coefficients = result["coefficients"]
    every_fold = coefficients["n_folds_nonzero"] == result["outer_splits"]
    nonzero = coefficients.loc[every_fold, "gene"].astype(str).tolist()
    fold_nonzero = [row["n_nonzero"] for row in result["elastic_net"]["folds"]]
    kept = {
        key: value
        for key, value in result.items()
        if key != "coefficients"
    }
    kept["elastic_net_mean_nonzero"] = float(np.mean(fold_nonzero))
    kept["n_genes_nonzero_every_fold"] = len(nonzero)
    kept["genes_nonzero_every_fold"] = nonzero
    return kept


def run_within_batch(expression: pd.DataFrame, groups: pd.Series, config: AnalysisConfig) -> dict:
    return within_batch_summary(nested_classification(expression, groups, config))


def analyse_external_cohort(
    train: pd.DataFrame,
    train_groups: pd.Series,
    discovery: pd.DataFrame,
    external: pd.DataFrame,
    external_groups: pd.Series,
    meta: dict,
    outdir,
    config: AnalysisConfig,
) -> dict:
    """QC, differential expression, transfer, and nested CV on one external cohort."""

    outdir = Path(outdir)
    figure_dir = outdir / "figures"
    table_dir = outdir / "tables"
    figure_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    quality = sample_quality(external, config.correlation_threshold)
    retained = quality["pass_correlation"].fillna(False)
    expression = external.loc[retained]
    groups = external_groups.loc[expression.index]
    counts = groups.value_counts()
    if counts.get(TNBC, 0) < config.min_samples_per_group or counts.get(
        "non-TNBC", 0
    ) < config.min_samples_per_group:
        raise RuntimeError(f"Too few external samples after filtering: {counts.to_dict()}")

    table = differential_expression(
        expression,
        groups,
        fdr_threshold=config.fdr_threshold,
        log2fc_threshold=config.log2fc_threshold,
    )
    transfer = transfer_elastic_net(train, train_groups, expression, groups, config)
    auc = plot_roc(
        transfer.pop("y_true"),
        transfer.pop("scores"),
        figure_dir / "tcga_external_roc.png",
        "TCGA-BRCA external test of the GSE76275 elastic net",
    )
    transfer["roc_auc_on_figure"] = auc
    within = run_within_batch(expression, groups, config)
    overlap = replication_table(
        discovery, table, config.fdr_threshold, config.log2fc_threshold
    )
    overlap.to_csv(table_dir / "gse76275_replication_in_tcga.csv", index=False)
    table.to_csv(table_dir / "tcga_differential_expression.csv", index=False)
    finite_fdr = table["fdr"].dropna()
    scale = value_scale_summary(expression)
    return {
        **meta,
        "n_samples_analysed": int(expression.shape[0]),
        "n_tnbc": int(counts[TNBC]),
        "n_non_tnbc": int(counts["non-TNBC"]),
        "n_samples_excluded": int((~retained).sum()),
        "correlation_threshold": config.correlation_threshold,
        "correlation_min": float(quality.loc[retained, "median_correlation"].min()),
        "correlation_median": float(quality.loc[retained, "median_correlation"].median()),
        "correlation_max": float(quality.loc[retained, "median_correlation"].max()),
        "n_genes": int(expression.shape[1]),
        "scale": scale,
        "differential_expression": {
            "n_genes_tested": int(table.shape[0]),
            "fdr_threshold": config.fdr_threshold,
            "log2fc_threshold": config.log2fc_threshold,
            "n_fdr_only": int((finite_fdr <= config.fdr_threshold).sum()),
            "n_significant": int(table["significant"].sum()),
            "n_higher_in_tnbc": int(((table["significant"]) & (table["higher_in"] == TNBC)).sum()),
            "n_higher_in_non_tnbc": int(
                ((table["significant"]) & (table["higher_in"] == "non-TNBC")).sum()
            ),
            "median_welch_df": float(table["welch_df"].median()),
        },
        "transfer": {key: value for key, value in transfer.items() if key != "genes"}
        | {"genes_head": transfer["genes"][:15], "n_features": transfer["n_features"]},
        "within_batch": within,
        "replication": replication_summary(overlap),
    }
