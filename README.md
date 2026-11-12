# SurvFitSummary

**Survival Extrapolation with Standard Parametric and Spline Models**

## Overview

`SurvFitSummary` is an R package for the survival analyses that inform the
inputs of cost-effectiveness models submitted to health technology
assessment bodies. It prepares time-to-event data from data frames or files
(CSV, Excel, SAS, SAS transport), summarises them with Kaplan-Meier
estimates, Cox proportional hazards models and smoothed empirical hazards,
assesses the proportional hazards assumption, and fits standard parametric
and Royston-Parmar spline models to each arm and to both arms with the arm
as a covariate, with tables of parameters, covariance matrices, Cholesky
factors, goodness of fit and predicted survival.

Version 2.0.0 is under development. The Excel and Word outputs will be
added in the next step; the functions of version 1.0.0 were replaced (see
`NEWS.md`).

## Installation

```r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/SurvFitSummary")
```

## Main Functions

| Function | Description |
|---|---|
| `processing_dataset()` | Renames the arm, time and event columns to `ARM`, `SURVTIME` and `EVENT`, keeps covariates and strata, checks the values, and sets the control arm as the first level of `ARM` |
| `generate_dummy_survival_data()` | Simulates exponential event and dropout times |
| `summarise_data()` | Numbers of patients, events and censored observations, and the median with its confidence interval, by arm |
| `km_table()` | Kaplan-Meier estimates with standard errors and confidence intervals, and numbers at risk, events and censored observations at regular times |
| `estimate_hazard_ratio()` | Hazard ratios from unadjusted, adjusted, unstratified and stratified Cox models |
| `test_proportional_hazards()` | Grambsch-Therneau test of proportional hazards |
| `estimate_empirical_hazard()` | Smoothed empirical hazard by arm (`bshazard`) |
| `fit_survival_models()` | Fits parametric models (12 distributions, the 7 standard models of NICE DSU TSD 14 by default) and spline models on the hazard, odds and normal scales with 1 to 3 knots (`flexsurv`) |
| `parameter_table()` | Parameter estimates on the natural and transformed scales |
| `covariance_table()` | Covariance matrices and Cholesky factors on the transformed and natural scales |
| `gof_table()` | Log-likelihood, AIC and BIC sorted by AIC |
| `predict_survival()` | Survival, hazard or cumulative hazard at given times, with the Kaplan-Meier estimate |
| `survival_summary()` | Median, mean and restricted mean survival |
| `spline_knot_table()` | Placement of the knots of the spline models |
| `plot_km()` | Kaplan-Meier curves with numbers at risk, each ending at the last observed time of its arm |
| `plot_empirical_hazard()` | Smoothed empirical hazard without fitted models |
| `plot_survival_fits()` | Kaplan-Meier curve with fitted survival curves extrapolated to a time horizon |
| `plot_hazard_fits()` | Smoothed empirical hazard with fitted hazards |
| `plot_distribution_diagnostics()` | Transformations of the Kaplan-Meier estimate and of the hazard that are linear for the Weibull, log-logistic, lognormal and Gompertz distributions |
| `plot_schoenfeld()` | Scaled Schoenfeld residuals with the log hazard ratio and a regression on time |
| `plot_qq_survival()` | Quantile-quantile plot of the survival times of two arms |

## Examples

```r
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
gof_table(fits)
parameter_table(fits)
covariance_table(fits, type = "cholesky", scale = "transformed")
predict_survival(fits, times = seq(0, 60, by = 6))
survival_summary(fits, time_horizon = 60)
spline_knot_table(fits)

plot_km(fit_data, conf_int = TRUE)
plot_survival_fits(fits, arms = "Control", time_horizon = 120)
plot_hazard_fits(fits, models = "Log-logistic")
plot_distribution_diagnostics(fit_data)
plot_schoenfeld(fit_data)
plot_qq_survival(fit_data)
```

## References

- Jackson C (2016). flexsurv: A platform for parametric survival modeling
  in R. *Journal of Statistical Software* 70(8):1-33.
  doi:10.18637/jss.v070.i08
- Grambsch PM, Therneau TM (1994). Proportional hazards tests and
  diagnostics based on weighted residuals. *Biometrika* 81(3):515-526.
  doi:10.1093/biomet/81.3.515
- Rebora P, Salim A, Reilly M (2014). bshazard: A flexible tool for
  nonparametric smoothing of the hazard function. *The R Journal*
  6(2):114-122. doi:10.32614/RJ-2014-028
- Royston P, Parmar MKB (2002). Flexible parametric proportional-hazards
  and proportional-odds models for censored survival data, with
  application to prognostic modelling and estimation of treatment effects.
  *Statistics in Medicine* 21(15):2175-2197. doi:10.1002/sim.1203

## License

GPL (>= 2)
