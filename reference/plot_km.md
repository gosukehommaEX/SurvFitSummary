# Kaplan-Meier Curves

Plots the Kaplan-Meier curve of each selected arm, with or without
pointwise confidence bands, and the numbers at risk below the plot. Each
curve ends at the largest observed time of its arm, and the axes start
at (0, 0).

## Usage

``` r
plot_km(
  data,
  arms = NULL,
  conf_int = FALSE,
  conf_level = 0.95,
  risk_table = TRUE,
  time_unit = c("months", "days", "weeks", "years"),
  y_label = "Survival probability",
  break_by = NULL,
  legend_position = c("bottom", "right", "top", "left", "none"),
  colours = NULL,
  linetypes = NULL,
  text_sizes = NULL
)
```

## Arguments

- data:

  Data frame returned by
  [`processing_dataset`](https://gosukehommaEX.github.io/SurvFitSummary/reference/processing_dataset.md).

- arms:

  `NULL` (the default, all arms in one plot) or the arms to plot, for
  example a single arm.

- conf_int:

  Logical; `TRUE` adds the pointwise confidence bands (log
  transformation, as in
  [`km_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/km_table.md)).
  The default is `FALSE`.

- conf_level:

  Confidence level, a single number in (0, 1). The default is 0.95.

- risk_table:

  Logical; `TRUE` (the default) adds the numbers at risk at the breaks
  of the time axis.

- time_unit:

  Time unit of the data, for the label of the time axis and the spacing
  of its breaks: `"months"` (the default), `"days"`, `"weeks"` or
  `"years"`.

- y_label:

  Label of the y axis, for example `"Overall survival"`. The default is
  `"Survival probability"`.

- break_by:

  `NULL` (the default, chosen from usual intervals of the time unit) or
  the spacing of the breaks of the time axis.

- legend_position:

  Position of the legend: `"bottom"` (the default), `"right"`, `"top"`,
  `"left"` or `"none"`.

- colours:

  `NULL` (the default colours of the arms, the same in all plots) or
  colours of the arms: a vector named by arm, for some or all of the
  plotted arms, or an unnamed vector with one colour per plotted arm, in
  the order of the arms. The confidence bands and the numbers at risk
  take the colours of their arms.

- linetypes:

  `NULL` (the default, solid lines) or line types of the arms, given as
  for `colours`: `"solid"`, `"dashed"`, `"dotted"`, `"dotdash"`,
  `"longdash"`, `"twodash"`, `"blank"` or a string of 2, 4, 6 or 8
  hexadecimal digits (for example `"44"`).

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

A ggplot object, or with `risk_table = TRUE` a patchwork object whose
first element is the plot and second the table.

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
plot_km(dat)

plot_km(dat, arms = "Treatment", conf_int = TRUE)

```
