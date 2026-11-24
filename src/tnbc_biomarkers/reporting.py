"""Plain-language results text built only from a completed analysis summary."""

from __future__ import annotations


def _fmt_p(value: float | None) -> str:
    if value is None:
        return "NA"
    if value == 0:
        return "0"
    if value < 0.001:
        return f"{value:.3e}"
    return f"{value:.3f}"


def _fmt(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}"


def _metric_sentence(label: str, block: dict) -> str:
    mean = block["mean"]
    std = block["std"]
    return (
        f"{label}: pooled ROC AUC {_fmt(block['pooled_roc_auc'])}; "
        f"mean fold ROC AUC {_fmt(mean['roc_auc'])} (SD {_fmt(std['roc_auc'])}); "
        f"mean balanced accuracy {_fmt(mean['balanced_accuracy'])} "
        f"(SD {_fmt(std['balanced_accuracy'])}); "
        f"mean accuracy {_fmt(mean['accuracy'])} (SD {_fmt(std['accuracy'])}); "
        f"mean F1 {_fmt(mean['f1'])} (SD {_fmt(std['f1'])})."
    )


def _hyperparameter_sentence(models: dict) -> str:
    folds = models["elastic_net"]["folds"]
    cs = sorted({row["C"] for row in folds})
    ratios = sorted({row["l1_ratio"] for row in folds})
    c_text = ", ".join(_fmt(value, 2) for value in cs)
    ratio_text = ", ".join(_fmt(value, 2) for value in ratios)
    return (
        f"Elastic-net hyperparameters selected by the inner loop were C {c_text} "
        f"and l1_ratio {ratio_text} across outer folds. "
        f"The grids were C {', '.join(_fmt(value, 2) for value in models['cs'])} "
        f"and l1_ratio {', '.join(_fmt(value, 2) for value in models['l1_ratios'])}."
    )


def _concordance_sentence(qc: dict) -> str:
    concordance = qc.get("receptor_concordance")
    if concordance is None:
        return "No sample had ER, PR, and HER2 all recorded."
    return (
        f"Among {qc['n_receptor_complete']} samples with ER, PR, and HER2 all recorded, "
        f"the fraction whose receptor pattern matched the submitted triple-negative field was "
        f"{_fmt(concordance)}."
    )


def _series_breakdown(mapping: dict) -> str:
    pieces = []
    for series, groups in mapping.items():
        detail = ", ".join(f"{count} {group}" for group, count in groups.items())
        pieces.append(f"{series} ({detail})")
    return "; ".join(pieces)


def _leading_gene_table(genes: list[dict]) -> str:
    lines = [
        "",
        "Genes with the smallest FDR:",
        "",
        "| Gene | log2 fold change | FDR | Higher in | Called | Probe |",
        "| --- | ---: | ---: | --- | --- | --- |",
    ]
    for gene in genes:
        called = "yes" if gene.get("significant") else "no"
        lines.append(
            f"| {gene['gene']} | {_fmt(gene['log2_fold_change'], 2)} | {_fmt_p(gene['fdr'])} | "
            f"{gene['higher_in']} | {called} | {gene.get('probe_id', '')} |"
        )
    return "\n".join(lines)


def _elastic_net_sentence(models: dict) -> str:
    genes = models.get("genes_nonzero_every_fold") or []
    sentence = (
        f"The elastic net kept a mean of {_fmt(models['elastic_net_mean_nonzero'], 1)} "
        f"non-zero coefficients per outer fold. "
        f"{models['n_genes_nonzero_every_fold']} genes were non-zero in every outer fold."
    )
    if genes:
        shown = ", ".join(genes[:15])
        if len(genes) > 15:
            shown = f"{shown} (15 of {len(genes)} listed)"
        sentence += f" Those genes, ordered by mean absolute coefficient, start: {shown}."
    return sentence


def _term_rows(terms: list[dict], limit: int = 5) -> str:
    if not terms:
        return "No significant terms were returned.\n"
    lines = ["| Source | Term | ID | p-value | Intersection | Term size |", "| --- | --- | --- | ---: | ---: | ---: |"]
    for term in terms[:limit]:
        name = str(term["name"]).replace("|", "/")
        lines.append(
            f"| {term['source']} | {name} | {term['native']} | {_fmt_p(term['p_value'])} | "
            f"{term['intersection_size']} | {term['term_size']} |"
        )
    return "\n".join(lines) + "\n"


def render_results(summary: dict) -> str:
    """Render the results section from ``summary`` produced by one run."""

    qc = summary["qc"]
    de = summary["differential_expression"]
    models = summary["classification"]
    enrichment = summary["enrichment"]
    scale = qc["scale"]
    lines = [
        f"The run used GEO {summary['accession']} ({', '.join(summary['subseries'])}) on {summary['platform']}.",
        (
            f"After probe filtering and gene collapse, {de['n_genes_tested']} genes were tested in "
            f"{qc['n_samples_analysed']} tumors ({qc['n_tnbc']} TNBC and {qc['n_non_tnbc']} non-TNBC). "
            f"{qc['n_samples_excluded']} sample(s) failed the median pairwise Pearson correlation "
            f"threshold of {qc['correlation_threshold']:.2f}. "
            f"The minimum, median, and maximum of the retained median correlations were "
            f"{_fmt(qc['correlation_min'])}, {_fmt(qc['correlation_median'])}, and {_fmt(qc['correlation_max'])}."
        ),
        (
            f"Expression values had minimum {_fmt(scale['min'], 2)}, median {_fmt(scale['median'], 2)}, "
            f"and maximum {_fmt(scale['max'], 2)}. "
            + (
                "They were analysed as the submitted log2 RMA intensities."
                if scale["assumed_log2"]
                else "The median was above 30, so the values were not treated as log2 intensities."
            )
        ),
        (
            f"At FDR <= {de['fdr_threshold']} and absolute log2 fold change >= {de['log2fc_threshold']}, "
            f"{de['n_significant']} genes were called ({de['n_higher_in_tnbc']} higher in TNBC and "
            f"{de['n_higher_in_non_tnbc']} higher in non-TNBC). "
            f"{de['n_fdr_only']} genes had FDR <= {de['fdr_threshold']} with no fold-change filter. "
            f"The median Welch–Satterthwaite degrees of freedom was {_fmt(de['median_welch_df'], 1)}."
        ),
        (
            f"The smallest FDR was {_fmt_p(de['top_gene']['fdr'])} for {de['top_gene']['gene']} "
            f"(log2 fold change {_fmt(de['top_gene']['log2_fold_change'], 2)}, "
            f"probe {de['top_gene'].get('probe_id', 'NA')})."
        ),
        _leading_gene_table(de["leading_genes"]),
        (
            "Classifier metrics use TNBC as the positive class. "
            "Feature selection and tuning were nested: each outer training fold ranked genes by the "
            f"absolute Welch t statistic, kept {models['n_features']} genes, and chose hyperparameters "
            f"by {models['inner_splits']}-fold stratified ROC AUC. "
            f"Outer evaluation used {models['outer_splits']} stratified folds "
            f"(random state {models['random_state']}). "
            "Logistic models used balanced class weights and a probability threshold of 0.5."
        ),
        _metric_sentence("Elastic net", models["elastic_net"]),
        _metric_sentence("L2 logistic regression", models["l2_logistic"]),
        _metric_sentence("Majority-class baseline", models["majority_baseline"]),
        _hyperparameter_sentence(models),
        _elastic_net_sentence(models),
        (
            "One permutation of the group labels, using the same nested elastic net, "
            f"had pooled ROC AUC {_fmt(models['permuted_label_pooled_roc_auc'])}."
        ),
        (
            f"PC1 and PC2 of the {qc['pca']['n_genes']} most variable genes explained "
            f"{qc['pca']['pc1_variance_fraction'] * 100:.1f}% and "
            f"{qc['pca']['pc2_variance_fraction'] * 100:.1f}% of the variance."
        ),
        _concordance_sentence(qc),
        (
            f"Pathway queries used g:Profiler {enrichment.get('version')} "
            f"(timestamp {enrichment.get('timestamp')}), sources {', '.join(enrichment['sources'])}, "
            f"domain scope {enrichment.get('domain_scope')}, and the g:SCS threshold. "
            f"Background fallback to the annotated genome: {enrichment.get('background_fallback')}. "
            f"The TNBC-higher query contained {enrichment['up']['n_genes_queried']} of "
            f"{enrichment['up']['n_genes_passing']} passing genes and returned "
            f"{enrichment['up']['n_terms']} terms. "
            f"The non-TNBC-higher query contained {enrichment['down']['n_genes_queried']} of "
            f"{enrichment['down']['n_genes_passing']} passing genes and returned "
            f"{enrichment['down']['n_terms']} terms."
        ),
        "",
        "Top terms for genes higher in TNBC:",
        "",
        _term_rows(enrichment["up"]["top"]),
        "Top terms for genes higher in non-TNBC:",
        "",
        _term_rows(enrichment["down"]["top"]),
        "Sample labels by GEO subseries: " + _series_breakdown(qc["group_by_series"]) + ".",
        f"Wall time for this run was {summary['elapsed_seconds']:.1f} seconds.",
    ]
    return "\n".join(lines).rstrip() + "\n"
