# Changelog

## 1.1.0 - 2026-11-24

- Add TCGA-BRCA from UCSC Xena (HiSeqV2, one RNA-seq matrix) as an external cohort. TNBC is defined by ER, PR, and HER2 immunohistochemistry.
- Train an elastic net on GSE76275 and test it on TCGA-BRCA after per-cohort z-scoring. Also run nested cross-validation inside TCGA-BRCA.
- Count how many GSE76275 calls keep their direction and false-discovery rate in TCGA-BRCA.
- State the GSE76275 nested-CV scores as a batch-confounded internal result.

## 1.0.0 - 2026-11-24

- Analyse public GEO series GSE76275, a SuperSeries of triple-negative and non-triple-negative breast tumors profiled on the Affymetrix HG-U133 Plus 2.0 array.
- Download series matrices and the GPL570 annotation by accession, with local caching.
- Apply sample-correlation quality control and a principal-component plot.
- Test gene-level differential expression with a log2 fold change, a Welch t-test, and Benjamini–Hochberg FDR.
- Write a volcano plot and a clustered heatmap.
- Fit an elastic-net logistic classifier with nested cross-validation, and compare it with an L2 logistic model and a majority-class baseline.
- Query pathway enrichment from g:Profiler for genes higher in each group.
- Add tests, a pinned environment, and GitHub Actions CI. CI uses a cached subset of the processed cohort rather than downloading GEO.
