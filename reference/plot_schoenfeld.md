# Scaled Schoenfeld Residuals of the Treatment Effect

For each arm other than the control arm, fits a Cox model to the data of
the two arms and plots the scaled Schoenfeld residuals of the treatment
effect
([`survival::cox.zph()`](https://rdrr.io/pkg/survival/man/cox.zph.html))
against time with:

- the estimated log hazard ratio of the Cox model (pink line);

- the regression of the residuals on time (black line), a natural cubic
  spline with 4 degrees of freedom as in
  [`survival::plot.cox.zph()`](https://rdrr.io/pkg/survival/man/plot.cox.zph.html),
  with its pointwise confidence interval (dotted lines);

- the log hazard ratio and the p-value of the Grambsch-Therneau test of
  [`test_proportional_hazards`](https://gosukehommaEX.github.io/SurvFitSummary/reference/test_proportional_hazards.md),
  in the title of the panel.

A flat regression line close to the log hazard ratio indicates that the
treatment effect does not depend on time, that is, proportional hazards.

## Usage

``` r
plot_schoenfeld(
  data,
  conf_level = 0.95,
  time_unit = c("months", "days", "weeks", "years"),
  legend_position = c("bottom", "right", "top", "left", "none"),
  colours = NULL,
  linetypes = NULL,
  text_sizes = NULL
)
```

## Arguments

- data:

  Data frame returned by
  [`processing_dataset`](https://gosukehommaEX.github.io/SurvFitSummary/reference/processing_dataset.md),
  with two or more treatment arms.

- conf_level:

  Confidence level of the interval of the regression line, a single
  number in (0, 1). The default is 0.95.

- time_unit:

  Time unit of the data: `"months"` (the default), `"days"`, `"weeks"`
  or `"years"`.

- legend_position:

  Position of the legend: `"bottom"` (the default), `"right"`, `"top"`,
  `"left"` or `"none"`.

- colours:

  `NULL` (the default) or colours of the elements of the plot:
  `"residuals"` (points, dark blue), `"log_hr"` (estimated log hazard
  ratio, pink), `"regression"` (regression on time, black) and
  `"confidence_interval"` (black). A named vector sets the colours of
  the elements it names; an unnamed vector gives the four colours in
  this order.

- linetypes:

  `NULL` (the default) or line types of `"log_hr"` (solid),
  `"regression"` (solid) and `"confidence_interval"` (dotted), given as
  for `colours` (see
  [`plot_km`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_km.md)
  for the line types).

- text_sizes:

  `NULL` (the default) or a named vector of font sizes in points, with
  some or all of the names `"title"` (titles of the plot and of its
  panels), `"axis_title"`, `"axis_text"`, `"legend"` and `"risk_table"`
  (numbers at risk and names of the arms below the plot; the title of
  the table is 1.125 times as large). The texts not named keep their
  default sizes: 14 points for the axis titles, 16.8 for the title of a
  plot, and 11.2 for the titles of panels, the axis text, the legend and
  the numbers at risk.

## Value

A ggplot object (one panel per comparison with the control arm, titled
with the comparison, the log hazard ratio and the p-value).

## Details

The scaled residuals are centred on the log hazard ratio (they estimate
the time-varying coefficient beta(t)), so the reference line is the
estimated log hazard ratio and not 0. The confidence interval of the
regression uses the variance returned by `cox.zph()`, as in
`plot.cox.zph()`, with the normal quantile of `conf_level`.

## References

Grambsch PM, Therneau TM (1994). Proportional hazards tests and
diagnostics based on weighted residuals. Biometrika 81(3):515-526.
[doi:10.1093/biomet/81.3.515](https://doi.org/10.1093/biomet/81.3.515)

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
plot_schoenfeld(dat)

```
