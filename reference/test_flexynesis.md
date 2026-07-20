# Predict with a trained flexynesis model on the test cohort

Exports the test cohort to the flat inference layout and runs flexynesis
in inference mode against the model + artifacts from
[`train_flexynesis`](https://fracalanca.github.io/Omics2Surv/reference/train_flexynesis.md).

## Usage

``` r
test_flexynesis(
  training_result,
  X_test_list,
  y_test,
  tmpdir = tempfile("flexynesis_infer_")
)
```

## Arguments

- training_result:

  List returned by
  [`train_flexynesis`](https://fracalanca.github.io/Omics2Surv/reference/train_flexynesis.md).

- X_test_list:

  Named list of test omics matrices (samples x features).

- y_test:

  Data frame of test survival data; row names are sample IDs. `clin.csv`
  is written from this so known labels are available, but it is not
  required for prediction itself.

- tmpdir:

  Working directory for the exported test data; fresh by default.

## Value

A list with `pred`: the test-set predictions, ordered to match
`rownames(y_test)`.

## See also

[`train_flexynesis`](https://fracalanca.github.io/Omics2Surv/reference/train_flexynesis.md),
[`run_flexynesis_inference`](https://fracalanca.github.io/Omics2Surv/reference/run_flexynesis_inference.md)
