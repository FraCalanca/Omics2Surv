# Multi-omics survival analysis

``` r

library(Omics2Surv)
```

## Introduction

This vignette walks through a complete **multi-omics** survival analysis
workflow, from data acquisition to patient stratification and functional
enrichment. As a working example we use the TCGA-KIPAN cohort (kidney
pan-cancer), retrieving RNA-seq, copy number variation (CNV) and DNA
methylation data together with the corresponding clinical information.

While three omics layers are downloaded and assembled into a single
`MultiAssayExperiment` (MAE), this vignette focuses on **combining
multiple omics** (RNA-seq and Methylation), using early, late, and joint
integration strategies. For a workflow that models each omic
independently, see the companion
[`vignette("single-omic-analysis")`](https://fracalanca.github.io/Omics2Surv/articles/single-omic-analysis.md).

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

data_expl= explore_pca(KIPANRNA_link, KIPANClin_link, metavars = c("histological_type", "gender", "radiation_therapy", "race", "ethnicity"))

data_expl$eigencorplot
```

## Data preselection

### Optional splitting

The cohort is split into a training and a test set. The same split is
shared across all omics layers so that the same patients end up in the
training set (and the same patients in the test set) for every modality.

``` r

KIPAN_split <- split_train_test(
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

KIPAN_MAE_test=build_MAE(omics_list = omics_test, KIPAN_split$test_clinical)
```

## Multi-omics modeling: RNA + Methylation

From this point on we focus on combining the RNA-seq and methylation
layers (`omic_list = c("mRNA", "meth")`) using several integration
strategies.

### Penalized regression - Early integration

In early integration, the features from the different omics are
concatenated before fitting a single penalized regression model.

``` r

#lasso
KIPAN_MAE_train <- train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "early",
  model_method = "cox_lasso",
  nfolds = 10,
  nlambda = 200
)

KIPAN_MAE_test <- predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "early",
  model_method = "cox_lasso"
)
```

``` r

#adaptive lasso
KIPAN_MAE_train <- train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "early",
  model_method = "cox_adaptive_lasso",
  nfolds = 10,
  nlambda = 200
)

KIPAN_MAE_test <- predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "early",
  model_method = "cox_adaptive_lasso"
)
```

``` r

#elastic net
KIPAN_MAE_train <- train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "early",
  model_method = "cox_elastic_net",
  nfolds= 10,
  nlambda = 200
)

KIPAN_MAE_test <- predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "early",
  model_method = "cox_elastic_net"
)
```

``` r

#aft pen
KIPAN_MAE_train <- train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "early",
  model_method = "aft_pen",
  nfolds = 10,
  nlambda= 200
)

KIPAN_MAE_test <- predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "early",
  model_method = "aft_pen"
)
```

### Penalized regression - Late integration

In late integration, a separate model is fitted on each omic and the
resulting predictions are then combined. We used the Robust Rank
Aggregation method.

``` r

#lasso
KIPAN_MAE_train <- train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "late",
  model_method = "cox_lasso"
)

KIPAN_MAE_test <- predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "late",
  model_method = "cox_lasso"
)
```

``` r

#adaptive lasso
KIPAN_MAE_train <- train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "late",
  model_method = "cox_adaptive_lasso"
)

KIPAN_MAE_test <- predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "late",
  model_method = "cox_adaptive_lasso"
)
```

``` r

#elastic net
KIPAN_MAE_train <- train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "late",
  model_method = "cox_elastic_net"
)

KIPAN_MAE_test <- predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "late",
  model_method = "cox_elastic_net"
)
```

``` r

#aft pen
KIPAN_MAE_train <- train_pen_reg_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "late",
  model_method = "aft_pen"
)

KIPAN_MAE_test <- predict_pen_reg_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  integration  = "late",
  model_method = "aft_pen"
)
```

### Joint integration

Joint integration methods model all omics simultaneously within a single
algorithm. Here we use cooperative regularized Cox regression
(`CoopCox`), a random-forest based approach for blockwise data
(`blockForest`), and a cooperative AFT model (`AFTCoop`).

``` r

KIPAN_MAE_train <- train_j_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  joint_method = "CoopCox",
  folds = 10,
  seed = 8290
)

KIPAN_MAE_test <- predict_j_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  joint_method = "CoopCox"
)
```

``` r

KIPAN_MAE_train <- train_j_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  joint_method = "blockForest",
  seed = 8290
)

KIPAN_MAE_test <- predict_j_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  joint_method = "blockForest"
)
```

``` r

KIPAN_MAE_train <- train_j_models(
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  joint_method = "AFTCoop"
)

KIPAN_MAE_test <- predict_j_models(
  mae_pred     = KIPAN_MAE_test,
  mae_train    = KIPAN_MAE_train,
  omic_list    = c("mRNA", "meth"),
  joint_method = "AFTCoop"
)
```

### Flexynesis

Finally, we train a deep-learning based multi-omics model using
[Flexynesis](https://github.com/BIMSBbioinfo/flexynesis). Flexynesis is
a Python package, so before running the chunks below you need to set up
a dedicated Python environment with Python \>= 3.11 and install
Flexynesis in it:

``` bash
# requires mamba (https://mamba.readthedocs.io/en/latest/installation.html)
mamba create --name flexynesisenv python==3.11
mamba activate flexynesisenv

# install the latest version from PyPI (https://pypi.org/project/flexynesis)
python -m pip install flexynesis --upgrade
```

``` r

KIPAN_MAE_train <- train_flexynesis_models(
  mae_train  = KIPAN_MAE_train,
  omic_list  = c("mRNA", "meth"),
  threads= 19 #change this to the number of threads available on your machine
)

KIPAN_MAE_test <- predict_flexynesis_models(
  mae_pred   = KIPAN_MAE_test,
  mae_train  = KIPAN_MAE_train,
  omic_list  = c("mRNA", "meth")
)
```

## C-index tables

[`benchmark_models()`](https://fracalanca.github.io/Omics2Surv/reference/benchmark_models.md)
collects the concordance index (C-index) achieved by every model trained
so far and returns a comparison table. Here we restrict the comparison
to the early, late, joint, and Flexynesis.

``` r

KIPAN_table= benchmark_models(KIPAN_MAE_train, KIPAN_MAE_test,
                 keep_integrations = c("early", "late", "joint", "flexynesis"))
```

## Patient stratification

Using the risk scores from the best performing multi-omics model
(`blockForest` in this example), patients can be stratified into risk
groups and the results visualized with Kaplan-Meier curves, UMAP
embeddings, and heatmaps.

``` r

KIPAN_multi_strat <- patient_stratification(KIPAN_MAE_train, KIPAN_MAE_test,
                                            model = "blockForest",
                                            method = "survcutpoint",
                                            plot_scores = TRUE)
plot_km(KIPAN_multi_strat)
plot_umap(KIPAN_MAE_test, KIPAN_multi_strat, model = "blockForest", mae_train = KIPAN_MAE_train, top_feat = 100)
```

``` r

plot_heatmap(KIPAN_MAE_test, KIPAN_multi_strat, model = "blockForest",
             top_feat = 50, mae_train = KIPAN_MAE_train)
```

Finally, the user defined top features by model importance can be tested
for functional enrichment.

``` r

KIPAN_block <- enrich_features(KIPAN_MAE_train, model = "blockForest",
                                  top_feat = 50,
                                  strip_prefix = "^(mRNA|meth)_",
                               pvalue_cutoff = 0.07,
  qvalue_cutoff = 0.2)
plot_enrichment(KIPAN_block, "GO")
```

## Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
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
#> [1] Omics2Surv_1.0.1
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.6         knitr_1.52        rlang_1.3.0       xfun_0.62        
#>  [5] otel_0.2.0        textshaping_1.0.5 jsonlite_2.0.0    htmltools_0.5.9  
#>  [9] ragg_1.5.2        sass_0.4.10       rmarkdown_2.32    grid_4.6.1       
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
    #> 6                glmnet     5.1
    #> 7                   M3C  1.34.0
    #> 8  MultiAssayExperiment  1.38.0
    #> 9              parallel   4.6.1
    #> 10             PCAtools  2.24.0
    #> 11               penAFT   0.3.2
    #> 12                readr   2.2.0
    #> 13            Rdimtools   1.1.4
    #> 14           reticulate  1.47.0
    #> 15     RobustRankAggreg   1.2.1
    #> 16            S4Vectors  0.50.3
    #> 17 SummarizedExperiment  1.42.0
    #> 18              survAUC   1.4-0
    #> 19             survival   3.8-6
    #> 20         TCGAbiolinks  2.40.0
    #> 21          blockForest   0.2.7
    #> 22      clusterProfiler  4.20.3
    #> 23           enrichplot  1.32.1
    #> 24             pheatmap  1.0.13
    #> 25           ReactomePA  1.56.0
    #> 26            survminer   0.5.2
    #> 27                 uwot   0.2.5
    #> 28                caret   7.0-1
