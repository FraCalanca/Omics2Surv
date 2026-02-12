
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Omics2Surv

<!-- badges: start -->

<!-- badges: end -->

**Omics2Surv** is an R package for **survival analysis with multi-omics
data**.

## 📦 Package overview

**Omics2Surv** is designed to:

- Handle heterogeneous multi-omics datasets

- Support multiple integration strategies:

\- None (single-omics modeling)

\- Early integration

\- Late integration

\- Joint / cooperative models

- Store all data in a standardized MultiAssayExperiment objects

- Enable fair model comparison using concordance index (C-index)

## Supported modeling strategies

1.  **Single-omics, early and late integration** models:

- Cox LASSO
- Cox Elastic Net
- Cox Adaptive LASSO
- Penalized AFT models

2.  **Joint integration** models:

- CoopCox (Cooperative Cox regression)
- AFTCoop (Cooperative AFT regression)
- blockForest (Random Forest)
- flexynesis (Python-based deep learning survival framework)

## Installation

You can install **Omics2Surv** from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("FraCalanca/Omics2Surv")
```

## 🏛 Funding

This work is supported by the PRIN 2022 PNRR P2022BLN38
project, Computational approaches for the integration of multi-omics
data funded by European Union - Next Generation EU, CUP B53D23027810001.
