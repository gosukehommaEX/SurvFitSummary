# Smoothed Empirical Hazard with Fitted Hazards

Plots the smoothed empirical hazard
([`estimate_empirical_hazard`](https://gosukehommaEX.github.io/SurvFitSummary/reference/estimate_empirical_hazard.md))
with the hazards of fitted parametric or spline models, extrapolated to
a time horizon. As in
[`plot_survival_fits`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_survival_fits.md),
either one arm with one or more models (the empirical hazard black, each
model with its colour and its AIC in the legend) or one model with one
or more arms (each arm with its colour, the empirical hazard solid and
the model dashed). The empirical hazard ends at the largest observed
time of each arm, and the axes start at (0, 0). For the empirical hazard
alone, use
[`plot_empirical_hazard`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_empirical_hazard.md).

## Usage

``` r
plot_hazard_fits(
  x,
  models = NULL,
  arms = NULL,
  fit_type = c("independent", "dependent"),
  conf_int = FALSE,
  time_horizon = NULL,
  y_max = NULL,
  y_label = "Hazard rate",
  break_by = NULL,
  n_points = 200,
  n_sim = 1000,
  seed = 1,
  legend_position = c("bottom", "right", "top", "left", "none"),
  colours = NULL,
  linetypes = NULL,
  text_sizes = NULL
)
```

## Arguments

- x:

  Object returned by
  [`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md).

- models:

  `NULL` (the default, all fitted models of `fit_type`) or model labels
  as in `x$models$model`, for example
  `c("Weibull", "Hazards, 2 knots")`, in the order of the legend.

- arms:

  `NULL` (the default, all arms) or the arms to plot.

- fit_type:

  `"independent"` (the default, the model fitted to each arm) or
  `"dependent"` (the model fitted to both arms with the arm as a
  covariate).

- conf_int:

  Logical; `TRUE` adds confidence bands (see Details). The default is
  `FALSE`.

- time_horizon:

  `NULL` (the default, the largest observed time of the plotted arms) or
  the largest time of the plot, to show the extrapolation.

- y_max:

  `NULL` (the default, see Details) or the upper limit of the y axis.

- y_label:

  Label of the y axis. The default is `"Hazard rate"`.

- break_by:

  `NULL` (the default) or the spacing of the breaks of the time axis.

- n_points:

  Number of times from 0 to the time horizon at which the models are
  evaluated. The default is 200.

- n_sim:

  Number of simulated parameter vectors for the confidence band of a
  model. The default is 1000.

- seed:

  Seed of the simulation. The default is 1.

- legend_position:

  Position of the legend: `"bottom"` (the default), `"right"`, `"top"`,
  `"left"` or `"none"`. Below or above the plot, the curves of one arm
  with several models are listed in two columns.

- colours:

  `NULL` (the default colours) or colours of the curves. With one arm,
  the curves are the smoothed empirical hazard (named
  `"Smoothed empirical hazard"`, black by default) and the models (named
  by their labels, as in `models`); with several arms, the curves are
  named by arm. A named vector sets the colours of the curves it names;
  an unnamed vector gives one colour per curve, in the order of the
  legend. With one arm, the confidence band of the model takes its
  colour; with several arms, the bands take the colours of the arms.

- linetypes:

  `NULL` (the default) or line types of the smoothed empirical hazard
  (`"Smoothed empirical hazard"`) and of the models (named by their
  labels), in both layouts, given as for `colours` (see
  [`plot_km`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_km.md)
  for the line types). By default the lines are solid, except that of
  the model with several arms (dashed).

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

A ggplot object.

## Details

The hazards of the models are evaluated at `n_points` times from
`time_horizon / n_points` to `time_horizon` (not at time 0, where some
hazards are infinite). By default the y axis goes up to 1.2 times the
largest of the empirical hazard (and its upper confidence limit with
`conf_int = TRUE`) and of the model hazards after 5\\ the time horizon,
so that the very large hazards of some models just after time 0 do not
set the scale. With `conf_int = TRUE` the confidence bands of the
empirical hazards are drawn, and, when one model is plotted, the
confidence band of the model simulated as in
[`predict_survival`](https://gosukehommaEX.github.io/SurvFitSummary/reference/predict_survival.md).

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(dat, distributions = c("exp", "weibull",
                                                   "llogis"),
                            spline_scales = NULL)
plot_hazard_fits(fits, arms = "Control", conf_int = TRUE)

```
