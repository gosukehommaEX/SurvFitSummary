# SurvFitSummary

**Survival Extrapolation with Standard Parametric and Spline Models**

## Overview

`SurvFitSummary` is an R package for the survival analyses that inform
the inputs of cost-effectiveness models submitted to health technology
assessment bodies. It prepares time-to-event data from data frames or
files (CSV, Excel, SAS, SAS transport), summarises them with
Kaplan-Meier estimates, Cox proportional hazards models and smoothed
empirical hazards, assesses the proportional hazards assumption, and
fits standard parametric and Royston-Parmar spline models to each arm
and to both arms with the arm as a covariate, with tables of parameters,
covariance matrices, Cholesky factors, goodness of fit and predicted
survival. The results are written to an Excel workbook and to a Word
report for health technology assessment (NICE DSU TSD 14 models and
spline-based extrapolation).

Version 2.0.0 replaces the functions of version 1.0.0 (see `NEWS.md`).
The documentation and the articles below are on the website
<https://gosukehommaEX.github.io/SurvFitSummary/>.

## Installation

``` r

# install.packages("devtools")
devtools::install_github("gosukehommaEX/SurvFitSummary",
                         build_vignettes = TRUE)
```

## Articles

| Vignette | Content |
|----|----|
| [`vignette("SurvFitSummary")`](https://gosukehommaEX.github.io/SurvFitSummary/articles/SurvFitSummary.md) | Getting started: data, Kaplan-Meier estimates, hazard ratio, diagnostics, model fitting, extrapolation, Excel workbook and Word report |
| [`vignette("models", package = "SurvFitSummary")`](https://gosukehommaEX.github.io/SurvFitSummary/articles/models.md) | The parametric and spline models: parameterisations, independent and dependent fits, knots, covariance matrices and Cholesky factors for probabilistic sensitivity analysis, predictions with confidence intervals, model comparison |
| [`vignette("validation", package = "SurvFitSummary")`](https://gosukehommaEX.github.io/SurvFitSummary/articles/validation.md) | The validation of the package and the reproduction of published results |

## Main Functions

| Function | Description |
|----|----|
| [`processing_dataset()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/processing_dataset.md) | Renames the arm, time and event columns to `ARM`, `SURVTIME` and `EVENT`, keeps covariates and strata, checks the values, and sets the control arm as the first level of `ARM` |
| [`generate_dummy_survival_data()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/generate_dummy_survival_data.md) | Simulates exponential event and dropout times |
| [`summarise_data()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/summarise_data.md) | Numbers of patients, events and censored observations, and the median with its confidence interval, by arm |
| [`km_table()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/km_table.md) | Kaplan-Meier estimates with standard errors and confidence intervals, and numbers at risk, events and censored observations at regular times |
| [`estimate_hazard_ratio()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/estimate_hazard_ratio.md) | Hazard ratios from unadjusted, adjusted, unstratified and stratified Cox models |
| [`test_proportional_hazards()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/test_proportional_hazards.md) | Grambsch-Therneau test of proportional hazards |
| [`estimate_empirical_hazard()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/estimate_empirical_hazard.md) | Smoothed empirical hazard by arm (`bshazard`) |
| [`fit_survival_models()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md) | Fits parametric models (12 distributions, the 7 standard models of NICE DSU TSD 14 by default) and spline models on the hazard, odds and normal scales with 1 to 3 knots (`flexsurv`) |
| [`summary()`](https://rdrr.io/r/base/summary.html), [`plot()`](https://rdrr.io/r/graphics/plot.default.html) | Methods for the fitted models: fit statistics and status of all the models, and the fitted survival curves or hazards |
| [`parameter_table()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/parameter_table.md) | Parameter estimates on the natural and transformed scales |
| [`covariance_table()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/covariance_table.md) | Covariance matrices and Cholesky factors on the transformed and natural scales |
| [`gof_table()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/gof_table.md) | Log-likelihood, AIC and BIC sorted by AIC |
| [`predict_survival()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/predict_survival.md) | Survival, hazard or cumulative hazard at given times, with the Kaplan-Meier estimate |
| [`survival_summary()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/survival_summary.md) | Median, mean and restricted mean survival |
| [`spline_knot_table()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/spline_knot_table.md) | Placement of the knots of the spline models |
| [`plot_km()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_km.md) | Kaplan-Meier curves with numbers at risk, each ending at the last observed time of its arm |
| [`plot_empirical_hazard()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_empirical_hazard.md) | Smoothed empirical hazard without fitted models |
| [`plot_survival_fits()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_survival_fits.md) | Kaplan-Meier curve with fitted survival curves extrapolated to a time horizon |
| [`plot_hazard_fits()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_hazard_fits.md) | Smoothed empirical hazard with fitted hazards |
| [`plot_distribution_diagnostics()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_distribution_diagnostics.md) | Transformations of the Kaplan-Meier estimate and of the hazard that are linear for the Weibull, log-logistic, lognormal and Gompertz distributions |
| [`plot_schoenfeld()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_schoenfeld.md) | Scaled Schoenfeld residuals with the log hazard ratio and a regression on time |
| [`plot_qq_survival()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_qq_survival.md) | Quantile-quantile plot of the survival times of two arms |
| [`write_survival_workbook()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/write_survival_workbook.md) | Excel workbook with the tables of the analysis, including the parameters with their covariance matrices |
| [`write_survival_report()`](https://gosukehommaEX.github.io/SurvFitSummary/reference/write_survival_report.md) | Word report with methods, diagnostics, model fits and extrapolations, with figures with and without confidence intervals |

The plotting functions place the legend with `legend_position` and set
the colours, line types and font sizes with `colours`, `linetypes` and
`text_sizes`.

## Examples

``` r

library(SurvFitSummary)
raw <- generate_dummy_survival_data(
  arm = c("Control", "Treatment"), n = c(100, 100),
  hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 123
)
dat <- processing_dataset(raw, column_arm = "ARM",
                          column_survtime = "SURVTIME",
                          column_cnsr = "CNSR", control_arm = "Control")
summarise_data(dat)
km_table(dat, interval = 6)
estimate_hazard_ratio(dat)
test_proportional_hazards(dat)

fit_data <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(fit_data)
summary(fits)
plot(fits, time_horizon = 120)
gof_table(fits)
parameter_table(fits)
covariance_table(fits, type = "cholesky", scale = "transformed")
predict_survival(fits, times = seq(0, 60, by = 6))
survival_summary(fits, time_horizon = 60)
spline_knot_table(fits)

plot_km(fit_data, conf_int = TRUE)
plot_km(fit_data, colours = c("black", "grey50"),
        linetypes = c("solid", "dashed"), text_sizes = c(axis_title = 16))
plot_survival_fits(fits, arms = "Control", time_horizon = 120)
plot_hazard_fits(fits, models = "Log-logistic")
plot_distribution_diagnostics(fit_data)
plot_schoenfeld(fit_data)
plot_qq_survival(fit_data)

write_survival_workbook(fits, "survival_results.xlsx")
write_survival_report(fits, "survival_report.docx",
                      outcome = "Overall survival")
```

## Validation

The unit tests compare every fitted model, table and output file with
reference values computed independently in Python (numpy and scipy,
without survival analysis libraries; `inst/validation/python`). The
scripts of `inst/validation` compare the package with closed forms,
[`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html),
[`rstpm2::stpm2()`](https://rdrr.io/pkg/rstpm2/man/gsm.html) and the
Python reference, check the Excel workbook and the Word report, and
estimate the coverage of the confidence intervals by simulation; the
scripts of `inst/reproduce` reproduce published results of Jackson
(2016) and Royston and Parmar (2002). Every comparison is also run with
the values of the package shifted, to show that it would detect an error
of that size. The verdicts are in `inst/validation/results` and
`inst/reproduce/results`, and are summarised in
[`vignette("validation", package = "SurvFitSummary")`](https://gosukehommaEX.github.io/SurvFitSummary/articles/validation.md).

## References

- Jackson C (2016). flexsurv: A platform for parametric survival
  modeling in R. *Journal of Statistical Software* 70(8):1-33.
  <doi:10.18637/jss.v070.i08>
- Grambsch PM, Therneau TM (1994). Proportional hazards tests and
  diagnostics based on weighted residuals. *Biometrika* 81(3):515-526.
  <doi:10.1093/biomet/81.3.515>
- Rebora P, Salim A, Reilly M (2014). bshazard: A flexible tool for
  nonparametric smoothing of the hazard function. *The R Journal*
  6(2):114-122. <doi:10.32614/RJ-2014-028>
- Royston P, Parmar MKB (2002). Flexible parametric proportional-hazards
  and proportional-odds models for censored survival data, with
  application to prognostic modelling and estimation of treatment
  effects. *Statistics in Medicine* 21(15):2175-2197.
  <doi:10.1002/sim.1203>

## License

GPL (\>= 2)
