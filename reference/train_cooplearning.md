# Train Cooperative Learning Model for Multi-Omics Survival Analysis

Trains a cooperative Lasso (CoopLasso) Cox regression model that
integrates multiple omics datasets for survival prediction. The method
encourages cooperation between different data modalities while
performing feature selection through L1 regularization.

## Usage

``` r
train_cooplearning(
  X_list,
  y_training,
  lambdaVector = seq(0.1, 0.9, 0.1),
  folds = 5,
  alpha = 1,
  maxit = 100,
  seed = 123
)
```

## Arguments

- X_list:

  A named list of 2 or 3 omics data matrices or data frames. Each
  element represents a different omics layer (e.g., transcriptomics,
  proteomics, methylation). All matrices must have:

  - Samples as rows with matching row names across all omics layers

  - Features as columns

- y_training:

  A data frame containing survival outcomes with:

  - Row names matching sample identifiers in `X_list`

  - Column `overall_survival`: Numeric survival time

  - Column `status`: Numeric event indicator (1 = event, 0 = censored)

- lambdaVector:

  Numeric vector of lambda values to evaluate during cross-validation.
  Default is `seq(0.1, 0.9, 0.1)`. Lambda controls the strength of
  regularization.

- folds:

  Integer specifying the number of cross-validation folds for lambda
  selection. Default is 5.

- alpha:

  Numeric value controlling the agreement between layers. `alpha = 1`
  for 2 layers, `alpha = 0.5` for 3 layers. Default is 1.

- maxit:

  Integer maximum number of iterations for the optimization algorithm.
  Default is 100.

- seed:

  Integer seed for reproducibility of train/validation split and
  cross-validation. Default is 123.

## Value

A list containing:

- model:

  The trained cooplasso model object containing coefficient vector `b`
  and other model parameters.

- selected_lambda:

  Numeric optimal lambda value selected through cross-validation.

- selected_features:

  Character vector of feature names with non-zero coefficients (selected
  features across all omics layers).

- pred:

  Numeric vector of risk scores for training samples. Higher values
  indicate higher risk (shorter predicted survival).

- cindex:

  Numeric concordance index (C-index) on training data, calculated with
  `reverse = TRUE` to account for risk score directionality.

- omics:

  Character vector of omics layer names from the input `X_list`.

## Details

The function performs the following workflow:

1.  Validates input structure (2 or 3 omics layers required)

2.  Aligns samples across all omics layers and survival data using
    common identifiers

3.  Splits data into training (90\\

4.  Performs k-fold cross-validation to select optimal lambda parameter

5.  Trains final cooperative Lasso model with selected lambda

6.  Computes risk scores and C-index on training data

7.  Identifies selected features with non-zero coefficients
