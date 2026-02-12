# Late Integration of Penalized Regression Models Across Multiple Omics Layers

## Usage

``` r
L_pen_regression_models(
  mae,
  omic_list,
  surv_list,
  method = "RRA",
  model = c("cox_lasso", "cox_adaptive_lasso", "cox_elastic_net", "aft_pen"),
  train_fraction = 0.7,
  seed = 123
)
```

## Arguments

- mae:

  A MultiAssayExperiment object containing multiple omics datasets and
  clinical information. Results will be stored in
  `metadata(mae)$analysis$late_integration`.

- omic_list:

  A list of omics layer names (or SummarizedExperiment objects)
  corresponding to assays in the MAE. Each will be modeled
  independently. Common samples across all layers are automatically
  identified and used.

- surv_list:

  A specification for extracting survival data from the MAE, typically
  containing column names for time and event variables. Passed to
  `prepare_surv_data`.

- method:

  Character string specifying the late integration method for combining
  predictions across omics layers. Must be one of:

  - `"average"`: Simple or weighted average of risk scores

  - `"RRA"`: Robust Rank Aggregation, combines rankings rather than raw
    scores (default)

  Default is `"average"`.

- model:

  Character string specifying the penalized regression model to fit on
  each omics layer. Must be one of:

  - `"cox_lasso"`: Cox regression with Lasso penalty

  - `"cox_adaptive_lasso"`: Cox regression with adaptive Lasso penalty

  - `"cox_elastic_net"`: Cox regression with elastic net penalty

  - `"aft_pen"`: Accelerated Failure Time model with elastic net or
    grouped lasso penalties

  Default is `"cox_lasso"`.

- train_fraction:

  Numeric value between 0 and 1 specifying the proportion of samples to
  use for training. The remainder is used for testing. Default is 0.7
  (70\\

  seedInteger seed for reproducibility of the train/test split. Default
  is 123. The same seed ensures identical splits across all omics
  layers.

The input MultiAssayExperiment object with updated metadata. Results are
stored at `metadata(mae)$analysis$late_integration[[model]]` containing:

- `train`: Late integration results on training data, including:

  - Ensemble predictions (combined risk scores)

  - Ensemble C-index

  - Individual omics-specific predictions and C-indices

  - Integration method used

- `test`: Late integration results on test data with the same structure

Additionally, individual omics results are stored at:
`metadata(mae)$analysis[[model]][[omic_name]]` Performs late integration
of survival predictions from penalized regression models trained
independently on each omics layer. Each omics dataset is modeled
separately using the same train/test split, and predictions are combined
using averaging or rank-based aggregation methods. The function
implements a late integration (late fusion) workflow:

1.  Identifies common samples across all omics layers

2.  Creates a single reference train/test split based on the first omics
    layer

3.  Applies the same split to all other omics layers (ensures
    consistency)

4.  Trains the specified penalized model independently on each omics
    layer

5.  Combines predictions from all layers using the specified integration
    method

6.  Evaluates ensemble performance on both training and test sets

7.  Stores individual and ensemble results in MAE metadata

Integration Methods

- average:

  Combines risk scores by averaging. Can be simple mean or weighted
  based on layer-specific performance. Assumes risk scores are on
  comparable scales.

- RRA:

  Robust Rank Aggregation converts predictions to ranks within each
  layer, then combines ranks. More robust to outliers and scale
  differences. Useful when absolute risk scores vary widely across
  layers.

Result OrganizationResults are stored hierarchically:

    metadata(mae)
      └─ analysis
         ├─ [model]                     # Individual omics results
         │  ├─ [omic_1]
         │  │  ├─ train
         │  │  └─ test
         │  └─ [omic_2]
         │     ├─ train
         │     └─ test
         └─ late_integration            # Ensemble results
            └─ [model]
               ├─ train (ensemble)
               └─ test (ensemble)

[`S_or_E_pen_regression_models`](https://fracalanca.github.io/Omics2Surv/reference/S_or_E_pen_regression_models.md)
for single omics penalized regression
[`late_integration`](https://fracalanca.github.io/Omics2Surv/reference/late_integration.md)
for the integration function
[`split_train_test`](https://fracalanca.github.io/Omics2Surv/reference/split_train_test.md)
for data splitting
[`prepare_design_matrix`](https://fracalanca.github.io/Omics2Surv/reference/prepare_design_matrix.md)
for omics data preparation
