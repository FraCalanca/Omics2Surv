# Test Cooperative Learning Model

Evaluates a trained cooperative learning model on test data by
generating risk predictions and calculating the concordance index.

## Usage

``` r
test_cooplearning(X_list_test, y_test, training_result)
```

## Arguments

- X_list_test:

  A named list of omics data matrices or data frames for testing. Must
  contain the same omics layers (in the same order and with matching
  feature names) as used in training. Each element should have:

  - Samples as rows with matching row names across all omics layers

  - Features as columns, matching those in the training data

- y_test:

  A data frame containing test survival outcomes with:

  - Row names matching sample identifiers in `X_list_test`

  - Column `overall_survival`: Numeric survival time

  - Column `status`: Numeric event indicator (1 = event, 0 = censored)

- training_result:

  A list object returned by
  [`train_cooplearning`](https://fracalanca.github.io/Omics2Surv/reference/train_cooplearning.md)
  containing the trained model and coefficients.

## Value

A list containing:

- pred:

  Numeric vector of risk scores for test samples. Higher values indicate
  higher risk (shorter predicted survival). Calculated as the linear
  combination of features weighted by learned coefficients.

- cindex:

  Numeric concordance index (C-index) on test data. Values range from 0
  to 1, where 0.5 indicates random predictions and 1.0 indicates perfect
  concordance. Calculated with `reverse = TRUE` to account for risk
  score directionality.

## Details

This function applies the trained cooperative learning model to new test
data:

1.  Validates that `X_list_test` is a list

2.  Aligns test samples across all omics layers using common identifiers

3.  Combines all omics layers into a single feature matrix

4.  Ensures feature names match the training model

5.  Computes risk scores as a weighted sum: `X_test %*% coefficients`

6.  Calculates C-index to evaluate prediction performance

The C-index measures the model's discriminative ability—the probability
that, for a randomly selected pair of patients, the model correctly
ranks which patient will experience the event first.

## Data Requirements

- Test omics layers must match training layers in number, order, and
  features

- Only samples present in all test omics layers are used

- Feature names must be consistent with training data

## See also

[`train_cooplearning`](https://fracalanca.github.io/Omics2Surv/reference/train_cooplearning.md)
for training the model
