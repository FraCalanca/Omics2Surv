# Writes one cohort to disk in the directory structure flexynesis expects.

Two layouts are supported:

- `"train_test"`: a `train/` and a `test/` sub-folder, each with one CSV
  per omics layer plus `clin.csv`. This is the layout for a normal
  training run (`--data_path`).

- `"flat"`: `clin.csv` and the omics CSVs written directly in `outdir`,
  with no sub-folders. This is the layout required by inference mode
  (`--data_path_test`).

## Usage

``` r
export_flexynesis_data(
  X_list,
  y,
  outdir,
  layout = c("train_test", "flat"),
  X_list2 = NULL,
  y2 = NULL
)
```

## Arguments

- X_list:

  Named list of omics matrices (samples x features) for this cohort.

- y:

  Data frame of clinical/survival outcome; row names are sample IDs.

- outdir:

  Output directory (created recursively if needed).

- layout:

  One of `"train_test"` or `"flat"`.

- X_list2, y2:

  Only for `layout = "train_test"`: the second cohort written into
  `test/` (the first goes into `train/`).

## Value

Invisibly the `outdir` path.

## Details

Omics matrices are transposed on write (features as rows, samples as
columns), matching the flexynesis input convention. Inputs are therefore
expected to be samples-by-features. Sample concordance across layers and
`y` is asserted first to prevent silent patient loss.
