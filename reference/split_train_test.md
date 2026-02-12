# Split Data into Training and Testing Sets

Randomly partitions a design matrix and corresponding survival data into
training and testing subsets based on a specified proportion.

## Usage

``` r
split_train_test(design_matrix, surv_data, seed = 4820, train_fraction = NULL)
```

## Arguments

- design_matrix:

  A matrix or data frame of predictor variables (omics data), where rows
  represent samples.

- surv_data:

  A matrix or data frame containing survival outcomes, with row names
  matching `design_matrix`.

- seed:

  An integer used to set the random seed for reproducibility. Default is
  4820.

- train_fraction:

  A numeric value between 0 and 1 indicating the proportion of samples
  to include in the training set.

## Value

A named list containing four elements:

- `train_omic`: Training subset of the design matrix.

- `test_omic`: Testing subset of the design matrix.

- `train_survival`: Training subset of the survival data.

- `test_survival`: Testing subset of the survival data.
