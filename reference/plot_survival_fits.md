# Kaplan-Meier Curves with Fitted Survival Curves

Plots the Kaplan-Meier curve with the survival curves of fitted
parametric or spline models, extrapolated to a time horizon, with the
numbers at risk below the plot. Two layouts are available:

- one arm with one or more models: the Kaplan-Meier curve is black and
  each model has its colour, with its AIC in the legend (for example
  `"Log-logistic (AIC = 1300.7)"`);

- one model with one or more arms: each arm has its colour, the
  Kaplan-Meier curve is solid and the model dashed.

The Kaplan-Meier curves end at the largest observed time of each arm,
and the axes start at (0, 0).

## Usage

``` r
plot_survival_fits(
  x,
  models = NULL,
  arms = NULL,
  fit_type = c("independent", "dependent"),
  conf_int = FALSE,
  time_horizon = NULL,
  risk_table = TRUE,
  y_label = "Survival probability",
  break_by = NULL,
  n_points = 200,
  n_sim = 1000,
  seed = 1
)
```

## Arguments

- x:

  Object returned by
  [`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md).

- models:

  `NULL` (the default, all fitted models of `fit_type`) or model labels
  as in `x$models$model`, for example
  `c("Weibull", "Hazards, 2 knots")`.

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

- risk_table:

  Logical; `TRUE` (the default) adds the numbers at risk at the breaks
  of the time axis.

- y_label:

  Label of the y axis, for example `"Overall survival"`. The default is
  `"Survival probability"`.

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

## Value

A ggplot object, or with `risk_table = TRUE` a patchwork object whose
first element is the plot and second the table.

## Details

With `conf_int = TRUE` the pointwise confidence bands of the
Kaplan-Meier curves are drawn, and, when one model is plotted, the
confidence band of the model simulated as in
[`predict_survival`](https://gosukehommaEX.github.io/SurvFitSummary/reference/predict_survival.md).
The plot is drawn from the fitted object, so plots with and without
confidence bands do not refit the models.

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
plot_survival_fits(fits, arms = "Control", time_horizon = 120)

plot_survival_fits(fits, models = "Weibull", fit_type = "dependent",
                   conf_int = TRUE)

```
