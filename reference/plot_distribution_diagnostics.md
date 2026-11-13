# Diagnostic Plots of the Parametric Distributions

Plots transformations of the Kaplan-Meier estimate (or of the smoothed
empirical hazard) of each arm that are straight lines when the data
follow a distribution:

- `"weibull"`: `log(-log(S(t)))` against `log(t)`, straight for the
  Weibull (slope = shape; slope 1 for the exponential). Parallel lines
  also indicate proportional hazards, so this is the log cumulative
  hazard plot of the proportional hazards assessment.

- `"loglogistic"`: `log((1 - S(t)) / S(t))` against `log(t)`, straight
  for the log-logistic.

- `"lognormal"`: `qnorm(1 - S(t))` against `log(t)`, straight for the
  lognormal.

- `"gompertz"`: `log(h(t))` of the smoothed empirical hazard against
  `t`, straight for the Gompertz.

The Kaplan-Meier estimate is used at the event times where it is
strictly between 0 and 1. With `regression = TRUE` the least-squares
line of each arm is added (dashed).

## Usage

``` r
plot_distribution_diagnostics(
  data,
  which = c("weibull", "loglogistic", "lognormal", "gompertz"),
  regression = TRUE,
  time_unit = c("months", "days", "weeks", "years")
)
```

## Arguments

- data:

  Data frame returned by
  [`processing_dataset`](https://gosukehommaEX.github.io/SurvFitSummary/reference/processing_dataset.md).

- which:

  Character vector with the plots, a subset of
  `c("weibull", "loglogistic", "lognormal", "gompertz")` (the default,
  all four).

- regression:

  Logical; `TRUE` (the default) adds the least-squares line of each arm.

- time_unit:

  Time unit of the data: `"months"` (the default), `"days"`, `"weeks"`
  or `"years"`.

## Value

A ggplot object for a single plot, otherwise a patchwork object with the
plots in the order of `which`, two per row, and one legend.

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
plot_distribution_diagnostics(dat)

plot_distribution_diagnostics(dat, which = "weibull")

```
