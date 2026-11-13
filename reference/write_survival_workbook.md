# Write the Results to an Excel Workbook

Writes the tables of the analysis to an Excel workbook (with the writexl
package), one table per sheet with bold column names, with unrounded
numbers (empty cells for missing values) and the print-ready names of
the models (for example `"Log-logistic"` or `"Hazards, 2 knots"`). The
sheets are:

- Settings:

  Time unit, confidence level, reporting times and the versions of R and
  of the packages.

- Data summary:

  [`summarise_data`](https://gosukehommaEX.github.io/SurvFitSummary/reference/summarise_data.md).

- KM table:

  [`km_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/km_table.md)
  at `interval`.

- Hazard ratios, PH test:

  `hazard_ratios` and
  [`test_proportional_hazards`](https://gosukehommaEX.github.io/SurvFitSummary/reference/test_proportional_hazards.md)
  (two or more arms).

- Model status:

  Each model with its status (`"ok"`, `"warning"` or `"error"`) and
  message.

- Fit statistics:

  [`gof_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/gof_table.md).

- Parameters:

  Parameters on the natural scale, with the rows of their
  variance-covariance matrix in the columns `cov_1`, `cov_2`, ...
  (column k is the covariance with the k-th parameter of the model).

- Parameters (transformed):

  The same on the scale on which the parameters are estimated (for
  probabilistic sensitivity analysis).

- Covariance, Cholesky:

  [`covariance_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/covariance_table.md)
  on both scales.

- Survival summary:

  [`survival_summary`](https://gosukehommaEX.github.io/SurvFitSummary/reference/survival_summary.md)
  with the restricted mean up to `time_horizon`.

- Landmark survival:

  [`predict_survival`](https://gosukehommaEX.github.io/SurvFitSummary/reference/predict_survival.md)
  at `0, interval, 2 * interval, ...` up to `time_horizon` and at the
  landmark times.

- Spline knots:

  [`spline_knot_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/spline_knot_table.md)
  (spline models only).

## Usage

``` r
write_survival_workbook(
  x,
  file,
  interval = NULL,
  landmarks = NULL,
  time_horizon = NULL,
  hazard_ratios = NULL,
  conf_int = TRUE,
  n_sim = 1000,
  seed = 1
)
```

## Arguments

- x:

  Object returned by
  [`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md).

- file:

  Path of the `.xlsx` file to write (replaced if it exists).

- interval:

  Spacing of the reporting times of the Kaplan-Meier table and of the
  predictions. The default `NULL` uses the spacing of the time axis of
  the plots of the data (for example 12 months).

- landmarks:

  Times at which the survival is predicted. The default `NULL` uses 6
  months and 1, 2, 5 and 10 years, in the time unit of the data.

- time_horizon:

  Time horizon of the predictions and of the restricted mean. The
  default `NULL` uses the last landmark time.

- hazard_ratios:

  Data frame returned by
  [`estimate_hazard_ratio`](https://gosukehommaEX.github.io/SurvFitSummary/reference/estimate_hazard_ratio.md)
  (for example with the covariates and strata chosen for the analysis).
  The default `NULL` uses the Cox model of the arm without adjustment.

- conf_int:

  Logical; `TRUE` (the default) adds the simulated confidence intervals
  of the predictions and of the median, mean and restricted mean, which
  takes longer.

- n_sim, seed:

  Number of simulations and seed of the confidence intervals (see
  [`predict_survival`](https://gosukehommaEX.github.io/SurvFitSummary/reference/predict_survival.md)).

## Value

The path `file`, invisibly.

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(dat, distributions = c("exp", "weibull"),
                            spline_scales = NULL)
path <- tempfile(fileext = ".xlsx")
write_survival_workbook(fits, path, conf_int = FALSE)
unlink(path)
```
