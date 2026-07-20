# Evaluate Cox GLMNET Model on Test Data

Predicts risk scores for a fitted Cox model using independent test data.

## Usage

``` r
test_cox_glmnet(x_test, training_result)
```

## Arguments

- x_test:

  A matrix or data frame of predictor variables for the test set.

- training_result:

  A list containing the trained model objects:

  - `cv_fit`: A fitted `cv.glmnet` object.

  - `selected_lambda`: The penalty parameter (\\\lambda\\) to use for
    prediction.

## Value

A list with predictions: `pred`: A matrix of predicted risk scores
(linear predictors) for the test set.

## See also

[`train_cox_glmnet_LASSO`](https://fracalanca.github.io/Omics2Surv/reference/train_cox_glmnet_LASSO.md),
[`train_cox_glmnet_ElasticNet`](https://fracalanca.github.io/Omics2Surv/reference/train_cox_glmnet_ElasticNet.md),
[`train_cox_glmnet_adaptiveLASSO`](https://fracalanca.github.io/Omics2Surv/reference/train_cox_glmnet_adaptiveLASSO.md)
