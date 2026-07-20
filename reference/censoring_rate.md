# Censoring Rate of a Survival Status Vector

Computes the proportion of censored observations, i.e. the fraction of
entries in `status` that are not equal to `event_code`.

## Usage

``` r
censoring_rate(status, event_code = 1)
```

## Arguments

- status:

  Numeric or integer vector of event indicators.

- event_code:

  Value in `status` that denotes an observed event. Default is `1`;
  every other value is treated as censored.

## Value

A single numeric value in `[0, 1]`: the censoring rate.
