# Train Cooperative Learning Model for Multi-Omics Survival Analysis

Trains a cooperative Lasso Cox regression model integrating 2 or 3 omics
layers. An internal validation subset is held out for lambda selection
via cross-validation. Can be called directly or via
[`train_j_models`](https://fracalanca.github.io/Omics2Surv/reference/train_j_models.md).

## Usage

``` r
train_cooplearning(
  X_list,
  y_training,
  val_fraction = 0.1,
  lambdaVector = seq(0.1, 0.9, 0.1),
  folds = 5,
  alpha = 1,
  maxit = 100,
  seed = 123
)
```

## Arguments

- X_list:

  Named list of 2 or 3 omics matrices (samples x features). Row names
  must be sample identifiers, consistent across all layers.

- y_training:

  Data frame with `overall_survival` and `status` columns. Row names
  must match sample identifiers in `X_list`.

- val_fraction:

  Numeric proportion of training samples to hold out as internal
  validation for lambda selection. Default is `0.1`.

- lambdaVector:

  Numeric vector of lambda values to evaluate. Default is
  `seq(0.1, 0.9, 0.1)`.

- folds:

  Integer number of cross-validation folds. Default is 5.

- alpha:

  Numeric cooperation parameter: 1 for 2 layers, 0.5 for 3. Default is
  1.

- maxit:

  Integer maximum iterations. Default is 1000.

- seed:

  Integer seed. Default is 123.

## Value

A list with:

- `model`: Trained cooplasso model.

- `selected_lambda`: Optimal lambda selected by CV.

- `selected_features`: Non-zero coefficient feature names.

- `pred`: Named numeric vector of risk scores for the full input set
  (train + internal validation reassembled in original order).

- `omics`: Names of omics layers used.

## References

Hahn G, Prokopenko D, Hecker J, Lutz SM, Mullin K, Sejour L, Hide W,
Vlachos I, DeSantis S, Tanzi RE, Lange C. Prediction of disease-free
survival for precision medicine using cooperative learning on multi-omic
data. *Briefings in Bioinformatics*, 25(4):bbae267, July 2024.
[doi:10.1093/bib/bbae267](https://doi.org/10.1093/bib/bbae267)

## See also

[`test_cooplearning`](https://fracalanca.github.io/Omics2Surv/reference/test_cooplearning.md),
[`train_j_models`](https://fracalanca.github.io/Omics2Surv/reference/train_j_models.md)
