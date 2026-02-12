# Test Block Random Forest Model

Evaluates a trained block random forest model on test data by generating
risk predictions and calculating the concordance index.

## Usage

``` r
test_blockForest(x_test, y_test, training_result, num_threads = 10)
```

## Arguments

- x_test:

  A matrix or data frame of test predictors. Must have the same features
  (columns) as `x_training` used in training, including matching column
  names and prefixes.

- y_test:

  A data frame containing test survival outcomes with columns:

  - `overall_survival`: Numeric survival time

  - `status`: Numeric event indicator (1 = event, 0 = censored)

- training_result:

  A list object returned by
  [`train_blockForest`](https://fracalanca.github.io/Omics2Surv/reference/train_blockForest.md)
  containing the trained model and associated parameters. Note: The
  function expects `training_result$model_BF`, but `train_blockForest`
  returns `training_result$model`. Ensure compatibility or adjust
  accordingly.

- num_threads:

  Integer specifying the number of parallel threads to use for
  prediction. Default is 10.

## Value

A list containing:

- pred:

  Numeric vector of linear predictors (risk scores) for test samples,
  extracted from the final column of the cumulative hazard function
  (CHF). Higher values indicate higher risk (shorter predicted
  survival).

- cindex:

  Numeric concordance index (C-index) on test data. Values range from 0
  to 1, where 0.5 indicates random predictions and 1.0 indicates perfect
  concordance. Calculated with `reverse = TRUE` to account for risk
  score directionality.

## Details

This function applies the trained block random forest model to new test
data:

1.  Prepares test data in the same format as training data

2.  Generates predictions using the trained model

3.  Extracts risk scores from the cumulative hazard function (CHF)

4.  Calculates the concordance index to evaluate prediction accuracy
