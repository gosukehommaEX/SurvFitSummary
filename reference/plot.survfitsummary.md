# Plot the Survival Models Fitted by fit_survival_models()

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method of
the object returned by
[`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md):
plots the fitted survival curves with
[`plot_survival_fits`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_survival_fits.md)
(`type = "survival"`, the default) or the fitted hazards with
[`plot_hazard_fits`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_hazard_fits.md)
(`type = "hazard"`). These functions plot either one arm with one or
more models or one model with one or more arms, so with neither `arms`
nor a single model given, the models of the control arm (the first arm)
are plotted.

## Usage

``` r
# S3 method for class 'survfitsummary'
plot(x, type = c("survival", "hazard"), models = NULL, arms = NULL, ...)
```

## Arguments

- x:

  Object returned by
  [`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md).

- type:

  `"survival"` (the default) or `"hazard"`.

- models:

  `NULL` (the default, all fitted models of `fit_type`) or model labels
  as in `x$models$model`.

- arms:

  `NULL` (the default: all arms with a single model, otherwise the
  control arm) or the arms to plot.

- ...:

  Other arguments of
  [`plot_survival_fits`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_survival_fits.md)
  or
  [`plot_hazard_fits`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_hazard_fits.md),
  for example `fit_type`, `conf_int`, `time_horizon`, `colours` or
  `legend_position`.

## Value

The plot returned by
[`plot_survival_fits`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_survival_fits.md)
or
[`plot_hazard_fits`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_hazard_fits.md)
(a ggplot or patchwork object), drawn when it is printed.

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
plot(fits, time_horizon = 120)

plot(fits, type = "hazard", models = "Weibull")

```
