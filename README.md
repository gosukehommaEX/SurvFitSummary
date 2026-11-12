# SurvFitSummary

**Survival Extrapolation with Standard Parametric and Spline Models**

## Overview

`SurvFitSummary` is an R package for the survival analyses that inform the
inputs of cost-effectiveness models submitted to health technology
assessment bodies. It prepares time-to-event data from data frames or files
(CSV, Excel, SAS, SAS transport), summarises them with Kaplan-Meier
estimates, Cox proportional hazards models and smoothed empirical hazards,
and assesses the proportional hazards assumption.

Version 2.0.0 is under development. The fitting of standard parametric and
Royston-Parmar spline models, the tables of parameters, covariance matrices
and Cholesky factors, the plots, and the Excel and Word outputs will be
added in the next steps; the functions of version 1.0.0 were replaced (see
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
```

## References

- Grambsch PM, Therneau TM (1994). Proportional hazards tests and
  diagnostics based on weighted residuals. *Biometrika* 81(3):515-526.
  doi:10.1093/biomet/81.3.515
- Rebora P, Salim A, Reilly M (2014). bshazard: A flexible tool for
  nonparametric smoothing of the hazard function. *The R Journal*
  6(2):114-122. doi:10.32614/RJ-2014-028

## License

GPL (>= 2)
