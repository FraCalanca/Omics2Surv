# Train and Test Joint Multi-Omics Models

Applies joint (integrative) multi-omics survival analysis approaches.
Trains a selected joint model on training data, evaluates it on test
data, and stores the results in the MultiAssayExperiment metadata for
downstream analysis.

## Usage

``` r
J_models(
  mae,
  split,
  joint_method = c("CoopCox", "AFTCoop", "blockForest", "flexynesis"),
  omics_list,
  seed = 123,
  ...
)
```

## Arguments

- mae:

  A MultiAssayExperiment object containing the multi-omics data and
  clinical information. Results will be stored in
  `metadata(mae)$analysis$cooperative`.

- split:

  A list object returned by a data splitting function, containing:

  - `X_list`: List of omics matrices (samples as rows, features as
    columns)

  - `train_survival`: Data frame with training survival outcomes
    (columns: `overall_survival`, `status`)

  - `test_survival`: Data frame with test survival outcomes

  Row names of survival data frames are used to identify train/test
  samples.

- joint_method:

  Character string specifying the joint modeling method to use. Must be
  one of:

  - `"CoopCox"`: Cooperative Cox regression (supports 2-3 omics layers)

  - `"AFTCoop"`: Cooperative Accelerated Failure Time model (requires
    exactly 2 omics layers)

  - `"blockForest"`: Block random forest (supports any number of omics
    layers)

  - `"flexynesis"`: Flexible multi-omics integration via neural networks

- omics_list:

  A list of SummarizedExperiment objects or omics data. Required
  specifically for the `"flexynesis"` method; may be NULL for other
  methods.

- seed:

  Integer seed for reproducibility. Default is 123. Passed to the
  training functions to ensure consistent results.

- ...:

  Additional arguments passed to the specific training function. These
  vary by method:

  - `CoopCox`: `lambdaVector`, `folds`, `alpha`, `maxit`

  - `AFTCoop`: `model`, `rho_values`, `ncore_max_rho`, `ncore_max_cv`

  - `blockForest`: `prefixes`, `threads`, `always.split.variables`,
    `n_top_genes`

  - `flexynesis`: method-specific parameters

## Value

The input MultiAssayExperiment object with updated metadata. Training
and test results are stored at:
`metadata(mae)$analysis$cooperative[[joint_method]]` which contains:

- `train`: List of training results including model object, predictions,
  C-index, and selected features

- `test`: List of test results including predictions and C-index

The exact structure of train/test results depends on the specific method
used.

## Details

The function performs a unified workflow for joint multi-omics modeling:

1.  Extracts training and test sample identifiers from the split object

2.  Subsets each omics layer into training and test sets

3.  Trains the selected joint model on training data

4.  Evaluates the model on held-out test data

5.  Stores both training and test results in the MAE metadata

Each joint method has specific requirements:

- CoopCox:

  Handles 2-3 omics layers; performs cooperative regularization with
  cross-validation for lambda selection

- AFTCoop:

  Requires exactly 2 omics layers; models survival time directly using
  accelerated failure time framework with cooperative penalties

- blockForest:

  Accepts any number of omics layers; organizes features into blocks for
  structured random forest modeling

- flexynesis:

  Flexible neural network-based integration; requires the `omics_list`
  parameter

## Method-Specific Requirements

- **AFTCoop**: Will stop with error if `length(X_train_list) != 2`

- **flexynesis**: Requires `omics_list` parameter to be specified

- All methods require consistent sample identifiers across omics layers

## Result Storage

Results are stored hierarchically in the MAE metadata:

    metadata(mae)
      └─ analysis
         └─ cooperative
            └─ [joint_method]
               ├─ train (training results)
               └─ test (test results)

## See also

[`train_cooplearning`](https://fracalanca.github.io/Omics2Surv/reference/train_cooplearning.md),
[`test_cooplearning`](https://fracalanca.github.io/Omics2Surv/reference/test_cooplearning.md)
for CoopCox
[`train_AFTCoop`](https://fracalanca.github.io/Omics2Surv/reference/train_AFTCoop.md),
[`test_AFTCoop`](https://fracalanca.github.io/Omics2Surv/reference/test_AFTCoop.md)
for AFTCoop
[`train_blockForest`](https://fracalanca.github.io/Omics2Surv/reference/train_blockForest.md),
[`test_blockForest`](https://fracalanca.github.io/Omics2Surv/reference/test_blockForest.md)
for blockForest
