# Prepare Test Omics Data as SummarizedExperiment

Aligns a test omics matrix to the feature set of the corresponding
training `SummarizedExperiment`, retaining only the features selected
during training and ordering rows consistently. The result is packaged
as a `SummarizedExperiment` inheriting the preprocessing metadata of the
training object.

## Usage

``` r
prepare_test_SE(test_data, train_se, on_constant = c("drop", "error", "keep"))
```

## Arguments

- test_data:

  Numeric matrix of omics data for the test set (features as rows,
  samples as columns). Row names must be set.

- train_se:

  A `SummarizedExperiment` object produced by one of the `prepare_*`
  functions (mRNA, miRNA, CNV, Methylation). The features retained in
  `assays(train_se)$counts` define the target feature set.

- on_constant:

  Character. What to do with test samples that are constant across the
  (remaining) training features: `"drop"` (default) removes them,
  `"error"` stops, `"keep"` retains them.

## Value

A `SummarizedExperiment` object with:

- `assays$counts`: Test matrix restricted to training features, in the
  same row order as `train_se`, with constant samples handled per
  `on_constant`.

- `metadata$preprocessing`: Inherited from `train_se`, with
  `split = "test"` and `n_features` updated.

## See also

[`prepare_RNA`](https://fracalanca.github.io/Omics2Surv/reference/prepare_RNA.md),
[`prepare_methyl`](https://fracalanca.github.io/Omics2Surv/reference/prepare_methyl.md),
[`prepare_CNV`](https://fracalanca.github.io/Omics2Surv/reference/prepare_CNV.md)
