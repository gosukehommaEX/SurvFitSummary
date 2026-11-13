# Quantile-Quantile Plot of Survival Times

For each arm other than the control arm, plots the quantiles of the
survival time of the arm against those of the control arm, both from the
Kaplan-Meier estimates
([`quantile()`](https://rdrr.io/r/stats/quantile.html) of
[`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)),
with the identity line (dashed) and the least-squares line through the
points (solid). Points on a straight line through the origin indicate
that the survival times of the arm are a multiple of those of the
control arm (an accelerated failure time effect); points on the identity
line indicate no difference. The axes go from 0 to the largest quantile,
with the same scale on both axes.

## Usage

``` r
plot_qq_survival(
  data,
  probs = seq(0.05, 0.95, by = 0.05),
  time_unit = c("months", "days", "weeks", "years"),
  legend_position = c("bottom", "right", "top", "left", "none")
)
```

## Arguments

- data:

  Data frame returned by
  [`processing_dataset`](https://gosukehommaEX.github.io/SurvFitSummary/reference/processing_dataset.md),
  with two or more treatment arms.

- probs:

  Probabilities of the quantiles, numbers in (0, 1). The default is
  `seq(0.05, 0.95, by = 0.05)`. Quantiles that are not reached by a
  Kaplan-Meier curve are left out.

- time_unit:

  Time unit of the data: `"months"` (the default), `"days"`, `"weeks"`
  or `"years"`.

- legend_position:

  Position of the legend: `"bottom"` (the default), `"right"`, `"top"`,
  `"left"` or `"none"`.

## Value

A ggplot object (one panel per comparison with the control arm).

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
plot_qq_survival(dat)

```
