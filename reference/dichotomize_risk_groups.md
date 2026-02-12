# Dichotomize Continuous Risk Scores into High and Low Risk Groups

Converts continuous risk predictions into binary risk groups (High/Low)
based on a cutoff threshold. The cutoff can be determined from training
data using median, or a custom value, and is then applied to test data
for risk stratification.

## Usage

``` r
dichotomize_risk_groups(
  pred_train,
  pred_test,
  method = "median",
  custom_cutoff = NULL
)
```

## Arguments

- pred_train:

  Numeric vector of risk predictions from the training set. Used to
  determine the cutoff threshold. NA values are removed when calculating
  median cutoffs.

- pred_test:

  Numeric vector of risk predictions from the test set (or any dataset
  to be stratified). These predictions will be classified as High or Low
  risk based on the cutoff derived from `pred_train`.

- method:

  Character string specifying the method to determine the cutoff
  threshold. Options are:

  - `"median"`: Use the median of training predictions

  - `"custom"`: Use a user-specified cutoff value (requires
    `custom_cutoff`)

  Default is `"median"`.

- custom_cutoff:

  Numeric value specifying a custom cutoff threshold. Only used when
  `method = "custom"`. Must be provided if custom method is selected.
  Default is NULL.

## Value

A list containing:

- groups:

  Factor vector with two levels (`"Low"`, `"High"`) indicating risk
  group assignment for each sample in `pred_test`. Samples with
  predictions above the cutoff are classified as "High" risk; those at
  or below are "Low" risk.

- cutoff:

  Numeric value of the threshold used for dichotomization. This allows
  users to see the exact cutoff applied and use it consistently across
  multiple datasets.

## Details

This function is commonly used in survival analysis to:

- Stratify patients into risk groups for Kaplan-Meier analysis

- Create interpretable clinical risk categories from continuous scores

- Enable comparison of survival curves between high and low risk groups

- Facilitate risk-based treatment decision making

Classification logic:

- `pred_test > cutoff` → "High" risk

- `pred_test <= cutoff` → "Low" risk

## Cutoff Methods

- Median:

  Robust to outliers; splits training samples into equal-sized groups.

- Custom:

  Allows domain-specific thresholds based on clinical relevance, prior
  studies, or optimization criteria.
