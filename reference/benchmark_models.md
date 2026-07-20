# Collect C-index metrics from a train/test MAE pair

Traverses `@metadata$analysis` of both MAE objects and assembles a tidy
data frame with one row per (model, integration, dataset) combination.
Handles all five integration branches (single-omics, early, late, joint,
flexynesis) and the extra `beta_rho_*` nesting of AFTCoop.

## Usage

``` r
benchmark_models(
  mae_train,
  mae_test,
  aftcoop_rho = NULL,
  keep_datasets = NULL,
  keep_integrations = NULL
)
```

## Arguments

- mae_train:

  MAE object for the training split.

- mae_test:

  MAE object for the test split.

- aftcoop_rho:

  Numeric scalar selecting which AFTCoop `beta_rho` variant to keep
  (e.g. `0.5`). Must match one of the keys stored as `beta_rho_<value>`
  in `metrics`. If `NULL` (default) all three variants are returned as
  separate rows.

- keep_datasets:

  Character vector of dataset values to retain (e.g. `"mRNA"`,
  `c("mRNA", "early_integrated")`, `"multiomics"`). `NULL` (default)
  keeps all datasets. Matched against the `dataset` column of the
  output, which mirrors the keys used in `@metadata$analysis`.

- keep_integrations:

  Character vector of integration types to retain (e.g. `"single"`,
  `c("single", "early")`). `NULL` (default) keeps all. Matched against
  the `integration` column (`"single"`, `"early"`, `"late"`, `"joint"`,
  `"flexynesis"`).

## Value

A `data.frame` with columns `model`, `integration`, `dataset`,
`harrell_train`, `harrell_test`, `uno_test`, `n_features`. `n_features`
is populated for single-omics and early-integration rows; it is `NA` for
late, joint and flexynesis rows. Other missing values are also `NA`.
