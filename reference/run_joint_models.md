# Run Joint Multi-Omics Survival Models

Performs end-to-end joint modeling workflow: builds a
MultiAssayExperiment from omics data, prepares design matrices, splits
data into train/test sets, and fits a specified joint integration model.

## Usage

``` r
run_joint_models(
  omics_list,
  clinical_data,
  surv_list,
  joint_method = c("CoopCox", "AFTCoop", "blockForest", "flexynesis"),
  train_fraction = 0.7,
  seed = 123,
  ...
)
```

## Arguments

- omics_list:

  Named list of SummarizedExperiment objects, each representing an omics
  layer (e.g., transcriptomics, proteomics).

- clinical_data:

  Matrix or data frame of clinical data (variables as rows, samples as
  columns).

- surv_list:

  Character vector of survival variable names in clinical data (e.g.,
  `c("overall_survival", "status")`).

- joint_method:

  Character string specifying the joint model: `"CoopCox"` (2-3 omics),
  `"AFTCoop"` (exactly 2 omics), `"blockForest"` (any number), or
  `"flexynesis"`.

- train_fraction:

  Numeric proportion (0-1) of samples for training. Default is 0.7.

- seed:

  Integer seed for reproducibility. Default is 123.

- ...:

  Additional arguments passed to the specific joint modeling function.

## Value

A MultiAssayExperiment object with results stored in
`metadata(mae)$analysis$cooperative[[joint_method]]` containing:

- `train`: Training results (model, predictions, C-index, selected
  features)

- `test`: Test results (predictions, C-index)
