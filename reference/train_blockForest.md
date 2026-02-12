# Train Block Random Forest for Survival Analysis

Trains a block random forest model for survival data using the
blockForest algorithm. Features are organized into blocks based on
specified prefixes (e.g., clinical, genomic, transcriptomic), allowing
for structured variable selection and importance ranking.

## Usage

``` r
train_blockForest(
  x_training,
  y_training,
  prefixes = NULL,
  threads = 10,
  always.split.variables = NULL,
  n_top_genes = 100,
  seed = 8240
)
```

## Arguments

- x_training:

  A matrix or data frame of training predictors.

- y_training:

  A data frame containing survival outcomes with columns:

  - `overall_survival`: Numeric survival time

  - `status`: Numeric event indicator (1 = event, 0 = censored)

- prefixes:

  Character vector of prefixes used to define feature blocks. Each
  prefix identifies a group of related variables (e.g.,
  `c("clinical_", "omic1_", "omic2_")`). This parameter is required and
  must not be NULL.

- threads:

  Integer specifying the number of parallel threads to use. Default is
  10.

- always.split.variables:

  Optional list of variable names that should be considered at every
  split. Useful for forcing important clinical variables to be
  evaluated. Must be a list; set to NULL to disable. Default is NULL.

- n_top_genes:

  Integer specifying the number of top important features to select
  based on variable importance. Default is 100.

- seed:

  Integer seed for reproducibility. Default is 8240.

## Value

A list containing:

- model:

  The trained blockForest model object.

- cindex:

  Numeric training concordance index (C-index), calculated as 1 -
  prediction error.

- pred:

  Numeric vector of linear predictors (prognostic scores) for training
  samples, extracted from the cumulative hazard function.

- selected_features:

  Character vector of the top N most important features based on
  impurity-based variable importance, where N = `n_top_genes`.

- extra:

  List containing additional information:

  - `variable_importance`: Named numeric vector of importance scores for
    all features

  - `blocks`: List of feature block indices as defined by prefixes

## Details

The function performs the following steps:

1.  Organizes features into blocks based on the specified prefixes

2.  Calculates block-specific mtry values using the square root rule:
    `sqrt(block_size)`

3.  Trains a survival random forest using the logrank splitting rule

4.  Computes variable importance using impurity reduction

5.  Extracts training predictions and performance metrics

6.  Selects top important features for downstream analysis
