# Split Multi-Omics Data into Training and Test Sets

Identifies the common sample set across all omics layers, performs a
single stratified random split on that common set, and returns the raw
matrices partitioned into training and test subsets. This ensures every
patient lands in the same split across all omics layers, which is a
prerequisite for multi-omics survival analysis.

## Usage

``` r
split_train_test(omics_list, clinical_data, train_fraction = 0.7, seed = 4820)
```

## Arguments

- omics_list:

  Named list of raw omics matrices (features as rows, samples as
  columns). Matrices do not need to share the same samples beforehand —
  the intersection is computed internally.

- clinical_data:

  Matrix or data frame of clinical data (variables as rows, samples as
  columns). Must contain `"overall_survival"` and `"status"` rows, and
  must cover all common samples.

- train_fraction:

  Numeric value in (0, 1) specifying the proportion of common samples to
  assign to training. Default is 0.7.

- seed:

  Integer seed for reproducibility. Default is 4820.

## Value

A named list with:

- `train_ids`: Character vector of training sample IDs.

- `test_ids`: Character vector of test sample IDs.

- `train_list`: Named list of raw matrices subset to training samples
  (same names as `omics_list`).

- `test_list`: Named list of raw matrices subset to test samples.

- `train_clinical`: Clinical data matrix subset to training samples
  (variables as rows, samples as columns).

- `test_clinical`: Clinical data matrix subset to test samples.

## Details

The split is performed once on the intersection of all sample sets and
then applied identically to every omic layer and to `clinical_data`.
This guarantees that:

- No patient appears in training for one omic and test for another.

- `build_MAE` on the resulting train/test lists will not drop any
  additional samples when computing the common set.

Samples present in some but not all omics are silently excluded from
both splits. A message reports the common sample count and how many
samples were excluded per layer.

## See also

[`prepare_RNA`](https://fracalanca.github.io/Omics2Surv/reference/prepare_RNA.md),
[`prepare_test_SE`](https://fracalanca.github.io/Omics2Surv/reference/prepare_test_SE.md),
[`build_MAE`](https://fracalanca.github.io/Omics2Surv/reference/build_MAE.md)
