# Smoothed Empirical Hazard

Plots the smoothed empirical hazard of each selected arm
([`estimate_empirical_hazard`](https://gosukehommaEX.github.io/SurvFitSummary/reference/estimate_empirical_hazard.md)),
without fitted models, with or without pointwise confidence bands. The
axes start at (0, 0).

## Usage

``` r
plot_empirical_hazard(
  data,
  arms = NULL,
  conf_int = FALSE,
  conf_level = 0.95,
  time_unit = c("months", "days", "weeks", "years"),
  y_label = "Hazard rate",
  break_by = NULL,
  legend_position = c("bottom", "right", "top", "left", "none")
)
```

## Arguments

- data:

  Data frame returned by
  [`processing_dataset`](https://gosukehommaEX.github.io/SurvFitSummary/reference/processing_dataset.md).

- arms:

  `NULL` (the default, all arms in one plot) or the arms to plot.

- conf_int:

  Logical; `TRUE` adds the pointwise confidence bands of
  [`bshazard::bshazard()`](https://rdrr.io/pkg/bshazard/man/bshazard.html).
  The default is `FALSE`.

- conf_level:

  Confidence level, a single number in (0, 1). The default is 0.95.

- time_unit:

  Time unit of the data: `"months"` (the default), `"days"`, `"weeks"`
  or `"years"`.

- y_label:

  Label of the y axis. The default is `"Hazard rate"`.

- break_by:

  `NULL` (the default) or the spacing of the breaks of the time axis.

- legend_position:

  Position of the legend: `"bottom"` (the default), `"right"`, `"top"`,
  `"left"` or `"none"`.

## Value

A ggplot object.

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
plot_empirical_hazard(dat, conf_int = TRUE)

```
