The run used GEO GSE76275 (GSE76124, GSE76274) on GPL570.
After probe filtering and gene collapse, 20848 genes were tested in 265 tumors (198 TNBC and 67 non-TNBC). 0 sample(s) failed the median pairwise Pearson correlation threshold of 0.80. The minimum, median, and maximum of the retained median correlations were 0.856, 0.911, and 0.930.
Expression values had minimum 3.43, median 7.73, and maximum 15.62. They were analysed as the submitted log2 RMA intensities.
At FDR <= 0.05 and absolute log2 fold change >= 1.0, 486 genes were called (188 higher in TNBC and 298 higher in non-TNBC). 12897 genes had FDR <= 0.05 with no fold-change filter. The median Welch–Satterthwaite degrees of freedom was 126.5.
The smallest FDR was 5.291e-71 for FZD9 (log2 fold change 1.55, probe 207639_at).

Genes with the smallest FDR:

| Gene | log2 fold change | FDR | Higher in | Called | Probe |
| --- | ---: | ---: | --- | --- | --- |
| FZD9 | 1.55 | 5.291e-71 | TNBC | yes | 207639_at |
| COL18A1-AS1 | -1.01 | 1.137e-68 | non-TNBC | yes | 1564017_at |
| SKA1 | 1.17 | 9.298e-64 | TNBC | yes | 217640_x_at |
| MIR924HG | 0.62 | 3.998e-54 | TNBC | no | 240365_at |
| GTSE1 | 1.09 | 1.684e-53 | TNBC | yes | 211040_x_at |
| ORC1 | 1.30 | 6.258e-52 | TNBC | yes | 205085_at |
| EN1 | 2.46 | 6.258e-52 | TNBC | yes | 220559_at |
| DBF4B | 1.06 | 1.518e-49 | TNBC | yes | 238508_at |
Classifier metrics use TNBC as the positive class. Feature selection and tuning were nested: each outer training fold ranked genes by the absolute Welch t statistic, kept 250 genes, and chose hyperparameters by 5-fold stratified ROC AUC. Outer evaluation used 5 stratified folds (random state 42). Logistic models used balanced class weights and a probability threshold of 0.5.
Elastic net: pooled ROC AUC 1.000; mean fold ROC AUC 1.000 (SD 0.000); mean balanced accuracy 1.000 (SD 0.000); mean accuracy 1.000 (SD 0.000); mean F1 1.000 (SD 0.000).
L2 logistic regression: pooled ROC AUC 1.000; mean fold ROC AUC 1.000 (SD 0.000); mean balanced accuracy 1.000 (SD 0.000); mean accuracy 1.000 (SD 0.000); mean F1 1.000 (SD 0.000).
Majority-class baseline: pooled ROC AUC 0.488; mean fold ROC AUC 0.500 (SD 0.000); mean balanced accuracy 0.500 (SD 0.000); mean accuracy 0.747 (SD 0.010); mean F1 0.855 (SD 0.007).
Elastic-net hyperparameters selected by the inner loop were C 0.01 and l1_ratio 0.15 across outer folds. The grids were C 0.01, 0.10, 1.00, 10.00, 100.00 and l1_ratio 0.15, 0.50, 0.85.
The elastic net kept a mean of 127.8 non-zero coefficients per outer fold. 95 genes were non-zero in every outer fold. Those genes, ordered by mean absolute coefficient, start: COL18A1-AS1, FZD9, C5orf49, LOC101930608, CACNG4, TMEM75, GRIK3, GTSE1, E2F1, HEXA-AS1, ATP11A-AS1, SLC25A21-AS1, MIR4435-2HG, DEFB132, DMPK (15 of 95 listed).
One permutation of the group labels, using the same nested elastic net, had pooled ROC AUC 0.476.
PC1 and PC2 of the 2000 most variable genes explained 19.5% and 11.3% of the variance.
Among 248 samples with ER, PR, and HER2 all recorded, the fraction whose receptor pattern matched the submitted triple-negative field was 1.000.
Pathway queries used g:Profiler e114_eg62_p19_9dc92ac8 (timestamp 2026-11-24T16:17:33.749516+00:00), sources GO:BP, KEGG, REAC, domain scope custom, and the g:SCS threshold. Background fallback to the annotated genome: False. The TNBC-higher query contained 188 of 188 passing genes and returned 104 terms. The non-TNBC-higher query contained 298 of 298 passing genes and returned 7 terms.

Top terms for genes higher in TNBC:

| Source | Term | ID | p-value | Intersection | Term size |
| --- | --- | --- | ---: | ---: | ---: |
| GO:BP | epithelium development | GO:0060429 | 1.698e-07 | 39 | 1137 |
| KEGG | Cell cycle | KEGG:04110 | 3.828e-07 | 13 | 156 |
| KEGG | Cadherin signaling | KEGG:04519 | 1.325e-06 | 18 | 359 |
| REAC | Formation of the cornified envelope | REAC:R-HSA-6809371 | 1.397e-06 | 11 | 98 |
| GO:BP | mitotic cell cycle | GO:0000278 | 6.367e-06 | 31 | 838 |

Top terms for genes higher in non-TNBC:

| Source | Term | ID | p-value | Intersection | Term size |
| --- | --- | --- | ---: | ---: | ---: |
| GO:BP | response to endogenous stimulus | GO:0009719 | 0.002 | 47 | 1386 |
| KEGG | AMPK signaling pathway | KEGG:04152 | 0.005 | 9 | 117 |
| KEGG | Hormone signaling | KEGG:04081 | 0.007 | 12 | 212 |
| GO:BP | cellular response to endogenous stimulus | GO:0071495 | 0.008 | 41 | 1167 |
| GO:BP | mammary gland alveolus development | GO:0060749 | 0.015 | 6 | 20 |

Sample labels by GEO subseries: GSE76124 (198 TNBC, 0 non-TNBC); GSE76274 (0 TNBC, 67 non-TNBC).
Wall time for this run was 28.7 seconds.
