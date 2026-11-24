"""Run the GSE76275 TNBC versus non-TNBC analysis."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from tnbc_biomarkers.classify import nested_classification
from tnbc_biomarkers.config import AnalysisConfig
from tnbc_biomarkers.differential import differential_expression, passing_genes
from tnbc_biomarkers.enrichment import gprofiler_enrich
from tnbc_biomarkers.geo import collapse_to_genes, fetch_cohort
from tnbc_biomarkers.plots import plot_heatmap, plot_pca, plot_volcano
from tnbc_biomarkers.qc import sample_quality, value_scale_summary
from tnbc_biomarkers.reporting import render_results

TNBC = "TNBC"
NON_TNBC = "non-TNBC"


def run(
    accession: str,
    cache_dir: Path,
    outdir: Path,
    config: AnalysisConfig | None = None,
) -> dict:
    """Fetch ``accession`` and write tables, figures, and a results summary."""

    started = time.perf_counter()
    config = config or AnalysisConfig(accession=accession)
    print(f"Fetching {accession}", flush=True)
    cohort = fetch_cohort(accession, cache_dir)
    summary = analyse_expression(
        expression_by_probe=cohort["expression"],
        samples=cohort["samples"],
        symbols=cohort["symbols"],
        outdir=outdir,
        config=config,
        enrich_fn=gprofiler_enrich,
    )
    summary.update(
        {
            "accession": cohort["accession"],
            "subseries": cohort["subseries"],
            "platform": cohort["platform"],
            "soft_url": cohort["soft_url"],
            "matrix_urls": cohort["matrix_urls"],
            "annotation_url": cohort["annotation_url"],
            "checksums": cohort["checksums"],
            "data_processing": cohort["data_processing"],
            "annotation_date": cohort["annotation_date"],
            "elapsed_seconds": time.perf_counter() - started,
        }
    )
    _write_summary(summary, outdir)
    print(render_results(summary), flush=True)
    return summary


def analyse_expression(
    expression_by_probe: pd.DataFrame,
    samples: pd.DataFrame,
    symbols: pd.Series,
    outdir: Path,
    config: AnalysisConfig,
    enrich_fn=gprofiler_enrich,
) -> dict:
    """Analyse an in-memory probe matrix. Used by the command line and by tests."""

    outdir.mkdir(parents=True, exist_ok=True)
    figure_dir = outdir / "figures"
    table_dir = outdir / "tables"
    figure_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    gene_by_sample, probe_for_gene, probe_counts = collapse_to_genes(expression_by_probe, symbols)
    expression = gene_by_sample.T
    expression.index.name = "sample_id"
    samples = samples.loc[expression.index].copy()
    labelled = samples["group"].isin([TNBC, NON_TNBC])
    n_unlabelled = int((~labelled).sum())
    samples = samples.loc[labelled]
    expression = expression.loc[samples.index]

    quality = sample_quality(expression, config.correlation_threshold)
    scale = value_scale_summary(expression)
    if not scale["assumed_log2"]:
        raise RuntimeError(
            "Expression median is above 30, so the matrix does not look like log2 RMA. "
            "Refusing to treat it as log2 fold change."
        )
    retained = quality["pass_correlation"].fillna(False)
    analysed_samples = samples.loc[retained]
    analysed = expression.loc[retained]
    counts = analysed_samples["group"].value_counts()
    too_small = counts.get(TNBC, 0) < config.min_samples_per_group or counts.get(
        NON_TNBC, 0
    ) < config.min_samples_per_group
    if too_small:
        raise RuntimeError(
            "Too few samples remain after labelling and correlation filtering: "
            f"{counts.to_dict()}"
        )

    print(
        f"Analysing {analysed.shape[0]} samples and {analysed.shape[1]} genes",
        flush=True,
    )
    table = differential_expression(
        analysed,
        analysed_samples["group"],
        probe_for_gene=probe_for_gene,
        fdr_threshold=config.fdr_threshold,
        log2fc_threshold=config.log2fc_threshold,
    )
    pca_stats = plot_pca(
        analysed,
        analysed_samples["group"],
        figure_dir / "pca_samples.png",
        config.pca_genes,
    )
    plot_volcano(
        table,
        figure_dir / "volcano.png",
        config.fdr_threshold,
        config.log2fc_threshold,
        config.volcano_labels,
    )
    heatmap_genes = _heatmap_genes(table, config.heatmap_genes)
    plot_heatmap(
        analysed,
        analysed_samples["group"],
        heatmap_genes,
        figure_dir / "heatmap_top_genes.png",
    )
    print("Fitting nested classifiers", flush=True)
    classification = nested_classification(analysed, analysed_samples["group"], config)
    coefficients = classification.pop("coefficients")
    permuted_auc = _permuted_auc(analysed, analysed_samples["group"], config)
    background = analysed.columns.astype(str).tolist()
    enrichment = _enrich(table, background, config, enrich_fn)

    qc_table = analysed_samples.join(quality, how="left")
    failed = samples.loc[~retained].join(quality, how="left")
    qc_out = pd.concat([qc_table, failed], axis=0)
    qc_out = qc_out.reset_index().rename(columns={"index": "sample_id"})
    if "sample_id" not in qc_out.columns:
        qc_out = qc_out.rename(columns={qc_out.columns[0]: "sample_id"})
    keep_columns = [
        "sample_id",
        "title",
        "series",
        "group",
        "tn_status_raw",
        "er",
        "pr",
        "her2",
        "receptor_group",
        "data_processing",
        "median_expression",
        "iqr",
        "median_correlation",
        "pass_correlation",
    ]
    qc_out[keep_columns].to_csv(table_dir / "sample_qc.csv", index=False)
    table.to_csv(table_dir / "differential_expression.csv", index=False)
    table.loc[table["significant"]].to_csv(table_dir / "significant_genes.csv", index=False)
    _write_model_tables(classification, coefficients, table_dir)
    enrichment["up"]["table"].to_csv(table_dir / "enrichment_higher_in_tnbc.csv", index=False)
    enrichment["down"]["table"].to_csv(table_dir / "enrichment_higher_in_non_tnbc.csv", index=False)

    both = samples["receptor_group"].notna()
    concordance = None
    if int(both.sum()):
        agreement = samples.loc[both, "receptor_group"] == samples.loc[both, "group"]
        concordance = float(agreement.mean())
    crosstab = pd.crosstab(analysed_samples["group"], analysed_samples["series"])
    finite_fdr = table["fdr"].dropna()
    top = table.iloc[0]
    leading_columns = [
        "gene",
        "probe_id",
        "log2_fold_change",
        "fdr",
        "higher_in",
        "significant",
    ]
    leading = table.head(8)[leading_columns].to_dict(orient="records")
    every_fold = coefficients["n_folds_nonzero"] == config.outer_splits
    nonzero = coefficients.loc[every_fold, "gene"].tolist()
    fold_nonzero = [row["n_nonzero"] for row in classification["elastic_net"]["folds"]]
    mean_nonzero = float(np.mean(fold_nonzero))
    summary = {
        "qc": {
            "n_samples_downloaded": int(len(expression)),
            "n_unlabelled": n_unlabelled,
            "n_samples_excluded": int((~retained).sum()),
            "n_samples_analysed": int(len(analysed)),
            "n_tnbc": int(counts[TNBC]),
            "n_non_tnbc": int(counts[NON_TNBC]),
            "correlation_threshold": config.correlation_threshold,
            "correlation_min": float(quality.loc[retained, "median_correlation"].min()),
            "correlation_median": float(quality.loc[retained, "median_correlation"].median()),
            "correlation_max": float(quality.loc[retained, "median_correlation"].max()),
            "scale": scale,
            "probe_counts": probe_counts,
            "pca": pca_stats,
            "receptor_concordance": concordance,
            "n_receptor_complete": int(both.sum()),
            "group_by_series": _crosstab_dict(crosstab),
            "raw_tn_status_counts": (
                samples["tn_status_raw"].fillna("missing").value_counts().to_dict()
            ),
        },
        "differential_expression": {
            "n_genes_tested": int(table.shape[0]),
            "fdr_threshold": config.fdr_threshold,
            "log2fc_threshold": config.log2fc_threshold,
            "n_fdr_only": int((finite_fdr <= config.fdr_threshold).sum()),
            "n_significant": int(table["significant"].sum()),
            "n_higher_in_tnbc": int(((table["significant"]) & (table["higher_in"] == TNBC)).sum()),
            "n_higher_in_non_tnbc": int(
                ((table["significant"]) & (table["higher_in"] == NON_TNBC)).sum()
            ),
            "median_welch_df": float(table["welch_df"].median()),
            "top_gene": {
                "gene": str(top["gene"]),
                "probe_id": str(top["probe_id"]),
                "log2_fold_change": float(top["log2_fold_change"]),
                "fdr": float(top["fdr"]),
                "higher_in": str(top["higher_in"]),
            },
            "leading_genes": _records(leading),
        },
        "classification": {
            **{key: value for key, value in classification.items() if key != "coefficients"},
            "elastic_net_mean_nonzero": mean_nonzero,
            "n_genes_nonzero_every_fold": len(nonzero),
            "genes_nonzero_every_fold": nonzero,
            "permuted_label_pooled_roc_auc": permuted_auc,
        },
        "enrichment": enrichment["summary"],
        "heatmap_genes": heatmap_genes,
    }
    summary["_expression"] = analysed
    summary["_groups"] = analysed_samples["group"]
    summary["_table"] = table
    return summary


def write_ci_subset(summary: dict, path: Path, n_per_group: int = 8, n_each_tail: int = 20) -> None:
    """Write a small real matrix for tests. This is not a second analysis."""

    expression: pd.DataFrame = summary["_expression"]
    groups: pd.Series = summary["_groups"]
    table: pd.DataFrame = summary["_table"]
    ranked = table.assign(abs_t=table["t_statistic"].abs()).sort_values(
        ["abs_t", "gene"], ascending=[False, True], kind="mergesort"
    )
    high = ranked["gene"].head(n_each_tail).tolist()
    low = ranked["gene"].tail(n_each_tail).tolist()
    genes = list(dict.fromkeys([*high, *low]))
    chosen_samples = []
    for group in (TNBC, NON_TNBC):
        ids = sorted(groups.index[groups == group])[:n_per_group]
        chosen_samples.extend(ids)
    subset = expression.loc[chosen_samples, genes].copy()
    subset.insert(0, "group", groups.loc[chosen_samples].to_numpy())
    subset.insert(0, "sample_id", chosen_samples)
    path.parent.mkdir(parents=True, exist_ok=True)
    subset.to_csv(path, index=False)
    provenance = {
        "source_accession": summary.get("accession", "GSE76275"),
        "platform": summary.get("platform"),
        "subseries": summary.get("subseries"),
        "n_samples": len(chosen_samples),
        "n_genes": len(genes),
        "sample_ids": chosen_samples,
        "genes": genes,
        "note": (
            "Cached subset of the gene-level matrix from the GSE76275 run. "
            "Continuous integration uses it to exercise the same functions. "
            "Reported cohort results come from the full run, not from this subset."
        ),
    }
    provenance_path = path.with_suffix(".provenance.json")
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")


def _permuted_auc(expression: pd.DataFrame, groups: pd.Series, config: AnalysisConfig) -> float:
    """Pooled ROC AUC after one permutation of the group labels.

    This is a negative control for the fitted elastic net, not a permutation p-value.
    """

    generator = np.random.default_rng(config.random_state + 1)
    shuffled = groups.copy()
    shuffled.iloc[:] = generator.permutation(groups.to_numpy())
    print("Fitting one label-permuted elastic net", flush=True)
    permuted = nested_classification(expression, shuffled, config)
    permuted.pop("coefficients")
    return float(permuted["elastic_net"]["pooled_roc_auc"])


def _heatmap_genes(table: pd.DataFrame, limit: int) -> list[str]:
    called = table.loc[table["significant"], "gene"].head(limit).tolist()
    if len(called) >= min(limit, 10):
        return called
    return table["gene"].head(limit).tolist()


def _enrich(table: pd.DataFrame, background: list[str], config: AnalysisConfig, enrich_fn) -> dict:
    directions = {
        "up": passing_genes(table, TNBC, config.enrichment_cap),
        "down": passing_genes(table, NON_TNBC, config.enrichment_cap),
    }
    n_passing = {
        "up": int(((table["significant"]) & (table["higher_in"] == TNBC)).sum()),
        "down": int(((table["significant"]) & (table["higher_in"] == NON_TNBC)).sum()),
    }
    stored = {}
    metas = []
    for key, genes in directions.items():
        print(f"Enrichment query for {key}: {len(genes)} genes", flush=True)
        terms, meta = enrich_fn(genes, background, config.enrichment_sources)
        terms.insert(0, "direction", TNBC if key == "up" else NON_TNBC)
        stored[key] = {"table": terms, "meta": meta, "genes": genes}
        metas.append(meta)
    version = next((meta.get("version") for meta in metas if meta.get("version")), None)
    timestamp = next((meta.get("timestamp") for meta in metas if meta.get("timestamp")), None)
    scope = next((meta.get("domain_scope") for meta in metas if meta.get("domain_scope")), None)
    fallback = any(bool(meta.get("background_fallback")) for meta in metas)
    summary = {
        "service": "g:Profiler",
        "version": version,
        "timestamp": timestamp,
        "domain_scope": scope,
        "background_fallback": fallback,
        "sources": list(config.enrichment_sources),
        "up": _direction_summary(stored["up"], n_passing["up"]),
        "down": _direction_summary(stored["down"], n_passing["down"]),
    }
    return {"summary": summary, **stored}


def _direction_summary(stored: dict, n_passing: int) -> dict:
    table = stored["table"]
    top = []
    if not table.empty:
        top = _records(table.head(8).to_dict(orient="records"))
    return {
        "n_genes_passing": n_passing,
        "n_genes_queried": len(stored["genes"]),
        "n_terms": int(table.shape[0]),
        "top": top,
    }


def _write_model_tables(classification: dict, coefficients: pd.DataFrame, table_dir: Path) -> None:
    frames = []
    for name in ("elastic_net", "l2_logistic", "majority_baseline"):
        frame = pd.DataFrame(classification[name]["folds"])
        frame.insert(0, "model", name)
        frames.append(frame)
    pd.concat(frames, ignore_index=True).to_csv(table_dir / "classifier_folds.csv", index=False)
    coefficients.to_csv(table_dir / "elastic_net_coefficients.csv", index=False)


def _write_summary(summary: dict, outdir: Path) -> None:
    public = {key: value for key, value in summary.items() if not key.startswith("_")}
    text = json.dumps(public, indent=2, default=_json_default)
    (outdir / "summary.json").write_text(text + "\n", encoding="utf-8")
    (outdir / "RESULTS.md").write_text(render_results(public), encoding="utf-8")


def _crosstab_dict(crosstab: pd.DataFrame) -> dict:
    converted = {}
    for key, row in crosstab.to_dict().items():
        converted[str(key)] = {str(inner): int(value) for inner, value in row.items()}
    return converted


def _records(rows: list[dict]) -> list[dict]:
    converted = []
    for row in rows:
        converted.append(
            {
                key: value if isinstance(value, str) else _json_default(value)
                for key, value in row.items()
            }
        )
    return converted


def _json_default(value):
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        if not np.isfinite(value):
            return None
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (pd.Timestamp,)):
        return str(value)
    return value
