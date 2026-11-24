# Changelog

## 1.0.0 - 2026-11-24

- Analyse public GEO series GSE76275, a SuperSeries of triple-negative and non-triple-negative breast tumors profiled on the Affymetrix HG-U133 Plus 2.0 array.
- Download series matrices and the GPL570 annotation by accession, with local caching.
- Apply sample-correlation quality control and a principal-component plot.
- Test gene-level differential expression with a log2 fold change, a Welch t-test, and Benjamini–Hochberg FDR.
- Write a volcano plot and a clustered heatmap.
- Fit an elastic-net logistic classifier with nested cross-validation, and compare it with an L2 logistic model and a majority-class baseline.
- Query pathway enrichment from g:Profiler for genes higher in each group.
- Add tests, a pinned environment, and GitHub Actions CI. CI uses a cached subset of the processed cohort rather than downloading GEO.
