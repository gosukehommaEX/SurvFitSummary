# Write the Survival Analysis Report to a Word Document

Writes a Word report on the extrapolation of survival data for health
technology assessment (the standard parametric models of NICE DSU TSD 14
and spline-based models), with print-ready names of the models:

1.  Overview: data, outcome, arms, time unit, reporting interval and the
    variables of the hazard ratios. Items not given are left as
    placeholders in square brackets.

2.  Methods: Kaplan-Meier analysis, hazard ratios, diagnostics (with a
    table of the shapes of the hazard functions), standard parametric
    models, spline-based extrapolation and software versions.

3.  Results: data summary, Kaplan-Meier curves (all arms and each arm)
    and tables, hazard ratios, smoothed empirical hazards, distribution
    diagnostics, proportional hazards assessment (Schoenfeld residuals,
    Grambsch-Therneau test, QQ plot), and for each arm and for both arms
    fitted together: the standard parametric models (overlaid plots, AIC
    and BIC, and hazard and survival plots of each model) and the spline
    models (placement of the knots, fit statistics, median, mean and
    predicted survival at the landmark times, and overlaid plots).

The figures are drawn as PNG files and embedded in the document; the
versions without and with confidence intervals are both included by
default, so that they do not have to be recomputed. The tables of fit
statistics and predicted survival of the spline models of both arms are
given for each arm, so that they fit the page.

## Usage

``` r
write_survival_report(
  x,
  file,
  data_label = NULL,
  outcome = NULL,
  author = NULL,
  version = NULL,
  interval = NULL,
  landmarks = NULL,
  time_horizon = NULL,
  hazard_ratios = NULL,
  conf_int = c(FALSE, TRUE),
  n_sim = 1000,
  seed = 1,
  template = NULL,
  dpi = 300
)
```

## Arguments

- x:

  Object returned by
  [`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md).

- file:

  Path of the `.docx` file to write (replaced if it exists).

- data_label, outcome, author, version:

  `NULL` or a string: the data (trial, data cut-off, population), the
  outcome (for example `"Overall survival"`), the author and the version
  of the report.

- interval:

  Spacing of the reporting times of the Kaplan-Meier tables. The default
  `NULL` uses the spacing of the time axis of the plots of the data (for
  example 12 months).

- landmarks:

  Times of the predicted survival of the spline models. The default
  `NULL` uses 6 months and 1, 2, 5 and 10 years, in the time unit of the
  data.

- time_horizon:

  Time to which the models are extrapolated in the plots. The default
  `NULL` uses the last landmark time.

- hazard_ratios:

  Data frame returned by
  [`estimate_hazard_ratio`](https://gosukehommaEX.github.io/SurvFitSummary/reference/estimate_hazard_ratio.md)
  (for example with the covariates and strata chosen for the analysis).
  The default `NULL` uses the Cox model of the arm without adjustment.

- conf_int:

  Logical vector: the versions of the figures without (`FALSE`) and with
  (`TRUE`) confidence intervals; both by default.

- n_sim, seed:

  Number of simulations and seed of the confidence intervals of the
  models (see
  [`predict_survival`](https://gosukehommaEX.github.io/SurvFitSummary/reference/predict_survival.md)).

- template:

  `NULL` or the path of a `.docx` file whose page set-up and styles are
  used (for example an empty document with the styles of the
  organisation). The default `NULL` uses the template of the package: A4
  page, Arial 10 pt, numbered headings in dark red, bold captions, and
  tables with a dark red header row (made by
  `inst/templates/make_report_template.py`). The styles `"Title"`,
  `"Subtitle"`, `"heading 1"` to `"heading 4"`, `"Caption"`,
  `"Table Caption"`, `"Figure"`, `"Table Text"` and the table style
  `"Report Table"` are used when the document defines them.

- dpi:

  Resolution of the figures. The default is 300.

## Value

The path `file`, invisibly.

## Examples

``` r
# \donttest{
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(dat, distributions = c("exp", "weibull"),
                            spline_scales = NULL,
                            fit_type = "independent")
path <- tempfile(fileext = ".docx")
write_survival_report(fits, path, outcome = "Overall survival",
                      conf_int = FALSE, dpi = 100)
unlink(path)
# }
```
