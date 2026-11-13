# Median, Mean and Restricted Mean Survival of the Fitted Models

Reports, for each fitted model and each arm it predicts, the median
survival time, the mean survival time and, when `time_horizon` is given,
the restricted mean survival time up to the time horizon, as computed by
`summary.flexsurvreg()`.

## Usage

``` r
survival_summary(
  x,
  time_horizon = NULL,
  conf_int = FALSE,
  n_sim = 1000,
  seed = 1
)
```

## Arguments

- x:

  Object returned by
  [`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md).

- time_horizon:

  `NULL` (the default) or a single positive time horizon for the
  restricted mean, in the time unit of the data.

- conf_int:

  Logical; `TRUE` to add simulated confidence limits. The default is
  `FALSE`.

- n_sim:

  Number of simulated parameter vectors for the confidence limits. The
  default is 1000.

- seed:

  Seed of the simulation, a single whole number. The default is 1.

## Value

A data frame with the columns `fit_type`, `arm` (the arm predicted),
`model`, `median`, `median_lcl`, `median_ucl`, `mean`, `mean_lcl`,
`mean_ucl`, `rmst`, `rmst_lcl` and `rmst_ucl` (`NA` without
`time_horizon`; the limits are `NA` without `conf_int`).

## Details

The mean is the integral of the fitted survival function from 0 to
infinity; it is infinite for heavy-tailed fits (for example a
log-logistic shape of 1 or less) and for a Gompertz model with a
negative shape, which leaves a proportion of patients who never have the
event. Values that are not finite or cannot be computed are reported as
`NA`. The confidence limits, when requested, are simulated as in
[`predict_survival`](https://gosukehommaEX.github.io/SurvFitSummary/reference/predict_survival.md).

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(dat, distributions = c("exp", "weibull"),
                            spline_scales = NULL)
survival_summary(fits, time_horizon = 60)
#>      fit_type       arm       model    median median_lcl median_ucl     mean
#> 1 independent   Control Exponential  8.477043         NA         NA 12.22979
#> 2 independent   Control     Weibull  9.967365         NA         NA 12.02017
#> 3 independent Treatment Exponential 13.948793         NA         NA 20.12385
#> 4 independent Treatment     Weibull 15.677837         NA         NA 19.19071
#> 5   dependent   Control Exponential  8.477043         NA         NA 12.22979
#> 6   dependent Treatment Exponential 13.948793         NA         NA 20.12385
#> 7   dependent   Control     Weibull  9.887639         NA         NA 12.00441
#> 8   dependent Treatment     Weibull 15.800605         NA         NA 19.18324
#>   mean_lcl mean_ucl     rmst rmst_lcl rmst_ucl
#> 1       NA       NA 12.13927       NA       NA
#> 2       NA       NA 12.01768       NA       NA
#> 3       NA       NA 19.10328       NA       NA
#> 4       NA       NA 18.98189       NA       NA
#> 5       NA       NA 12.13927       NA       NA
#> 6       NA       NA 19.10328       NA       NA
#> 7       NA       NA 12.00129       NA       NA
#> 8       NA       NA 18.99793       NA       NA
```
