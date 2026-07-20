# Train Block Random Forest for Multi-Omics Survival Analysis

Trains a block random forest survival model. Feature blocks are derived
automatically from the names of `X_list`, prefixing each omic's column
names. Can be called directly or via
[`train_j_models`](https://fracalanca.github.io/Omics2Surv/reference/train_j_models.md).

## Usage

``` r
train_blockForest(
  X_list,
  y_training,
  threads = 10,
  always.split.variables = NULL,
  n_top_features = 100,
  seed = 123,
  ...
)
```

## Arguments

- X_list:

  Named list of omics matrices (samples x features). Names are used to
  derive block prefixes automatically (e.g. a layer named `"mRNA"` maps
  to columns starting with `"mRNA_"`). Column names in each matrix must
  already carry the prefix `"<omic_name>_"` (set automatically by
  [`train_j_models`](https://fracalanca.github.io/Omics2Surv/reference/train_j_models.md);
  must be set manually when calling directly).

- y_training:

  Data frame with `overall_survival` and `status`.

- threads:

  Integer number of parallel threads. Default is 10.

- always.split.variables:

  Optional list of variable names considered at every split. Default is
  `NULL`.

- n_top_features:

  Integer number of top important features to return. Default is 100.

- seed:

  Integer seed. Default is 123.

- ...:

  Additional arguments passed to
  [`blockForest::blockForest`](https://bips-hb.github.io/blockForest/reference/blockForest.html).

## Value

A list with:

- `model`: Trained blockForest model.

- `pred`: Numeric vector of risk scores for training samples.

- `selected_features`: Top `n_top_features` by impurity.

- `extra`: List with `variable_importance` and `blocks`.

## Details

Column names are prefixed internally before concatenation so that block
membership can be recovered from column names. If calling directly
(without
[`train_j_models`](https://fracalanca.github.io/Omics2Surv/reference/train_j_models.md)),
ensure column names in each matrix of `X_list` do not already carry the
prefix — they will be added here.

## See also

[`test_blockForest`](https://fracalanca.github.io/Omics2Surv/reference/test_blockForest.md),
[`train_j_models`](https://fracalanca.github.io/Omics2Surv/reference/train_j_models.md)
