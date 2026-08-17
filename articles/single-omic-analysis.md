# Single-omics survival analysis

``` r

library(Omics2Surv)
```

## Introduction

This vignette walks through a complete **single-omics** survival
analysis workflow, from data acquisition to patient stratification and
functional enrichment. As a working example we use the TCGA-KIPAN cohort
(Pan-Kidney Cohort), retrieving RNA-seq, copy number variance (CNV) and
DNA methylation data together with the corresponding clinical
information.

Although three omics layers are downloaded and assembled into a single
`MultiAssayExperiment` (MAE), in this vignette each omic is modeled
**independently**, one at a time. For a workflow that combines multiple
omics into a single predictive model, see the companion
[`vignette("multi-omic-analysis")`](https://fracalanca.github.io/Omics2Surv/articles/multi-omic-analysis.md).

## Data acquisition & exploration

We start by retrieving the omics and clinical data from LinkedOmics for
the TCGA-KIPAN cohort.

``` r

KIPANRNA_link = retrieve_linkedomics_data("TCGA-KIPAN","RNAseq (HiSeq, Gene level)")

KIPANCN_link = retrieve_linkedomics_data("TCGA-KIPAN", "SCNV (Gene level, Thresholded)")

KIPAN_Meth_link= retrieve_linkedomics_data("TCGA-KIPAN", "Methylation (Gene level, HM450K)")

KIPANClin_link= retrieve_linkedomics_data("TCGA-KIPAN", "Clinical")
```

### Exploration

Before any modeling step, it is good practice to explore the data
structure.
[`explore_pca()`](https://fracalanca.github.io/Omics2Surv/reference/explore_pca.md)
performs a principal component analysis on the RNA-seq data and relates
the resulting components to a set of clinical/metadata variables.

``` r

data_expl= explore_pca(KIPANRNA_link, KIPANClin_link, 
  metavars = c("histological_type", "gender", "radiation_therapy", "race", "ethnicity"))

data_expl$eigencorplot
```

## Data preselection

### Optional splitting

The cohort is split into a training and a test set. The same split is
shared across all omics layers so that the same patients end up in the
training set (and the same patients in the test set) for every layer.

``` r

KIPAN_split = split_train_test(
  omics_list     = list(mRNA = KIPANRNA_link,
                        CNV  = KIPANCN_link, methyl = KIPAN_Meth_link),
  clinical_data  = KIPANClin_link
)
```

### Variable screening

Each omic is preprocessed and screened on the training set only, to
avoid information leakage from the test set. Here RNA-seq features are
selected in an unsupervised fashion based on highly variable genes, CNV
features are selected based on sparsity, and methylation features are
selected based on interquartile range (IQR).

``` r

KIPAN_RNA_train=prepare_RNA(KIPAN_split$train_list$mRNA,   KIPAN_split$train_clinical,method = "unsupervised", unsupervised_method= "A", hvg_threshold = 5)
KIPAN_CNV_train=prepare_CNV(KIPAN_split$train_list$CNV,   KIPAN_split$train_clinical,method = "unsupervised", unsupervised_method = "sparsity", top_n_features = 800)
KIPAN_Meth_train=prepare_methyl(KIPAN_split$train_list$methyl,   KIPAN_split$train_clinical, method="unsupervised", unsupervised_method = "IQR", hvg_threshold = 5)
```

The same features selected on the training set are then applied to the
test set with
[`prepare_test_SE()`](https://fracalanca.github.io/Omics2Surv/reference/prepare_test_SE.md),
so that train and test share an identical feature space.

``` r

KIPAN_RNA_test= prepare_test_SE(KIPAN_split$test_list$mRNA,   KIPAN_RNA_train)
KIPAN_Meth_test= prepare_test_SE(KIPAN_split$test_list$methyl,  KIPAN_Meth_train)
KIPAN_CNV_test= prepare_test_SE(KIPAN_split$test_list$CNV,   KIPAN_CNV_train)
```

``` r

omics_train = list(mRNA=KIPAN_RNA_train, CNV=KIPAN_CNV_train, meth= KIPAN_Meth_train)

omics_test = list(mRNA= KIPAN_RNA_test, CNV= KIPAN_CNV_test, meth= KIPAN_Meth_test)
```

### MAE building with 3 omics

The preprocessed omics and the clinical data are assembled into a
`MultiAssayExperiment` object, separately for the training and the test
set.

``` r

KIPAN_MAE_train = build_MAE(omics_list = omics_train, KIPAN_split$train_clinical)

KIPAN_MAE_test = build_MAE(omics_list = omics_test, KIPAN_split$test_clinical)
```

Combine train and test clinical data to inspect the overall cohort, and
check the censoring rate of the survival status.

``` r

KIPAN_train_clin = as.data.frame(KIPAN_MAE_train@colData)
KIPAN_test_clin = as.data.frame(KIPAN_MAE_test@colData)
KIPAN_all_clin = as.data.frame(rbind(KIPAN_train_clin, KIPAN_test_clin))

censoring_rate(KIPAN_all_clin$status)

display_counts(KIPAN_all_clin,"status")
```

## Single-omics modeling

Even though the MAE contains three omics layers, the models below are
trained on **one omic at a time** (`integration = "none"`), so that each
modality’s predictive performance can be assessed independently.

### Penalized regression models

#### RNA

Four penalized regression models are trained on the RNA-seq data: a Cox
LASSO, a Cox Adaptive LASSO, a Cox Elastic Net, and a penalized
Accelerated Failure Time (AFT) model.

``` r

#cox_lasso
KIPAN_MAE_train = train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA"),
  integration  = "none",
  model_method = "cox_lasso",
  nlambda= 200,
  nfolds= 10,
  seed = 8290
)

KIPAN_MAE_test = predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA"),
  integration  = "none",
  model_method = "cox_lasso"
)
```

``` r

#adaptive lasso
KIPAN_MAE_train = train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA"),
  integration  = "none",
  model_method = "cox_adaptive_lasso",
  nlambda= 100,
  nfolds= 10,
  seed = 8290
)

KIPAN_MAE_test = predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA"),
  integration  = "none",
  model_method = "cox_adaptive_lasso"
)
```

``` r

#elastic net
KIPAN_MAE_train = train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA"),
  integration  = "none",
  model_method = "cox_elastic_net",
  nlambda= 200,
  nfolds= 10,
  seed = 8290
)

KIPAN_MAE_test = predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA"),
  integration  = "none",
  model_method = "cox_elastic_net"
)
```

``` r

#penAFT
KIPAN_MAE_train = train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA"),
  integration  = "none",
  model_method = "aft_pen",
  nlambda= 200,
  nfolds= 10,
  seed = 8290
)

KIPAN_MAE_test = predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA"),
  integration  = "none",
  model_method = "aft_pen"
)
```

#### Methylation

The same four models are trained on the methylation data.

``` r

#cox_lasso 
KIPAN_MAE_train = train_pen_reg_models( mae_train = KIPAN_MAE_train, omic_list = c("meth"), integration = "none", model_method = "cox_lasso", nlambda= 200, nfolds= 10, seed = 8290 )

KIPAN_MAE_test = predict_pen_reg_models( mae_pred = KIPAN_MAE_test, mae_train = KIPAN_MAE_train, omic_list = c("meth"), integration = "none", model_method = "cox_lasso" )
```

``` r

#adaptive lasso
KIPAN_MAE_train = train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("meth"),
  integration  = "none",
  model_method = "cox_adaptive_lasso",
  nlambda= 200,
  nfolds= 10,
  seed = 8290
)

KIPAN_MAE_test = predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("meth"),
  integration  = "none",
  model_method = "cox_adaptive_lasso"
)
```

``` r

#elastic net
KIPAN_MAE_train = train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("meth"),
  integration  = "none",
  model_method = "cox_elastic_net",
  nlambda= 200,
  nfolds= 10,
  seed = 8290
)

KIPAN_MAE_test = predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("meth"),
  integration  = "none",
  model_method = "cox_elastic_net"
)
```

``` r

#penAFT
KIPAN_MAE_train = train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("meth"),
  integration  = "none",
  model_method = "aft_pen",
  nlambda= 200,
  nfolds= 10,
  seed = 8290
)

KIPAN_MAE_test = predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("meth"),
  integration  = "none",
  model_method = "aft_pen"
)
```

#### CNV

Finally, the same four models are trained on the CNV data.

``` r

#cox_lasso
KIPAN_MAE_train = train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("CNV"),
  integration  = "none",
  model_method = "cox_lasso",
  nlambda= 200,
  nfolds= 10,
  seed = 8290
)

KIPAN_MAE_test = predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("CNV"),
  integration  = "none",
  model_method = "cox_lasso"
)
```

``` r

#adaptive lasso
KIPAN_MAE_train = train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("CNV"),
  integration  = "none",
  model_method = "cox_adaptive_lasso",
  nlambda= 200,
  nfolds= 10,
  seed = 8290
)

KIPAN_MAE_test = predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("CNV"),
  integration  = "none",
  model_method = "cox_adaptive_lasso"
)
```

``` r

#elastic net
KIPAN_MAE_train = train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("CNV"),
  integration  = "none",
  model_method = "cox_elastic_net",
  nlambda= 200,
  nfolds= 10,
  seed = 8290
)

KIPAN_MAE_test = predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("CNV"),
  integration  = "none",
  model_method = "cox_elastic_net"
)
```

``` r

KIPAN_MAE_train = train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("CNV"),
  integration  = "none",
  model_method = "aft_pen",
  nlambda= 200,
  nfolds= 10,
  seed = 8290
)

KIPAN_MAE_test = predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("CNV"),
  integration  = "none",
  model_method = "aft_pen"
)
```

## C-index tables

[`benchmark_models()`](https://fracalanca.github.io/Omics2Surv/reference/benchmark_models.md)
collects the concordance index (C-index) achieved by every model trained
so far and returns a comparison table. Here we restrict the comparison
to the three single-omic datasets.

``` r

KIPAN_single= benchmark_models(KIPAN_MAE_train, KIPAN_MAE_test, keep_datasets = c("mRNA","meth","CNV"))
```

## Patient stratification

Using the risk scores from the best performing single-omic model
(`aft_pen` on RNA in this example), patients can be stratified into risk
groups and the results visualized with Kaplan-Meier curves, UMAP
embeddings, and heatmaps.

``` r

KIPAN_strat = patient_stratification(KIPAN_MAE_train, KIPAN_MAE_test,
                               model = "aft_pen",                                modality = "mRNA",
                                 method = "survcutpoint",
                                 plot_scores = TRUE)
plot_km(KIPAN_strat)
plot_umap(KIPAN_MAE_test, KIPAN_strat, model = "aft_pen",
             modality = "mRNA", mae_train = KIPAN_MAE_train, n_neighbors = 10)
```

``` r

plot_heatmap(KIPAN_MAE_test, KIPAN_strat, model = "aft_pen",
             modality = "mRNA", mae_train = KIPAN_MAE_train)
```

Finally, the genes selected by the model can be tested for functional
enrichment.

``` r

KIPAN_RNA_en = enrich_features(KIPAN_MAE_train, model = "aft_pen", modality = "mRNA")
plot_enrichment(KIPAN_RNA_en, "GO")
```

## Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] Omics2Surv_1.0.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.6         knitr_1.51        rlang_1.3.0       xfun_0.60        
#>  [5] otel_0.2.0        textshaping_1.0.5 jsonlite_2.0.0    htmltools_0.5.9  
#>  [9] ragg_1.5.2        sass_0.4.10       rmarkdown_2.31    grid_4.6.1       
#> [13] evaluate_1.0.5    jquerylib_0.1.4   fastmap_1.2.0     yaml_2.3.12      
#> [17] lifecycle_1.0.5   compiler_4.6.1    fs_2.1.0          htmlwidgets_1.6.4
#> [21] lattice_0.22-9    systemfonts_1.3.2 digest_0.6.39     R6_2.6.1         
#> [25] splines_4.6.1     bslib_0.12.0      Matrix_1.7-5      tools_4.6.1      
#> [29] survival_3.8-6    pkgdown_2.2.1     cachem_1.1.0      desc_1.4.3
```

    #>                 package version
    #> 1               methods   4.6.1
    #> 2          BiocParallel  1.46.0
    #> 3            doParallel  1.0.17
    #> 4                 dplyr   1.2.1
    #> 5               ggplot2   4.0.3
    #> 6                glmnet     5.0
    #> 7                   M3C  1.34.0
    #> 8  MultiAssayExperiment  1.38.0
    #> 9              parallel   4.6.1
    #> 10             PCAtools  2.24.0
    #> 11               penAFT   0.3.2
    #> 12                readr   2.2.0
    #> 13            Rdimtools   1.1.4
    #> 14           reticulate  1.46.0
    #> 15     RobustRankAggreg   1.2.1
    #> 16            S4Vectors  0.50.1
    #> 17 SummarizedExperiment  1.42.0
    #> 18              survAUC   1.4-0
    #> 19             survival   3.8-6
    #> 20         TCGAbiolinks  2.40.0
    #> 21          blockForest   0.2.7
    #> 22      clusterProfiler  4.20.0
    #> 23           enrichplot  1.32.0
    #> 24         org.Hs.eg.db  3.23.1
    #> 25             pheatmap  1.0.13
    #> 26           ReactomePA  1.56.0
    #> 27            survminer   0.5.2
    #> 28                 uwot   0.2.4
