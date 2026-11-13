# Goodness-of-Fit Statistics of the Fitted Models

Reports the log-likelihood, the Akaike information criterion (AIC) and
the Bayesian information criterion (BIC) of each fitted model, sorted by
AIC in ascending order within each fit (each arm of the independent
fits, then the dependent fit).

## Usage

``` r
gof_table(x)
```

## Arguments

- x:

  Object returned by
  [`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md).

## Value

A data frame with the columns `fit_type`, `arm`, `model`, `n`, `events`,
`npars`, `loglik`, `AIC` and `BIC`. Models that failed are left out.

## Details

`AIC = -2 * loglik + 2 * npars` and
`BIC = -2 * loglik + log(n) * npars`, where `npars` is the number of
estimated parameters (the fixed shape of the Rayleigh model is not
counted) and `n` is the number of patients in the data fitted, censored
or not (the default of
[`flexsurv::BIC.flexsurvreg()`](http://chjackson.github.io/flexsurv-dev/reference/BIC.flexsurvreg.md)).
The log-likelihood is that of flexsurv, on the time scale (it includes
the log density of the event times).

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
gof_table(fits)
#>      fit_type       arm        model   n events npars     loglik      AIC
#> 1 independent   Control      Weibull 300    271     2  -931.5666 1867.133
#> 2 independent   Control Log-logistic 300    271     2  -936.5067 1877.013
#> 3 independent   Control  Exponential 300    271     1  -949.5500 1901.100
#> 4 independent Treatment      Weibull 300    240     2  -948.2785 1900.557
#> 5 independent Treatment Log-logistic 300    240     2  -950.0985 1904.197
#> 6 independent Treatment  Exponential 300    240     1  -960.4574 1922.915
#> 7   dependent Both arms      Weibull 600    511     3 -1879.9570 3765.914
#> 8   dependent Both arms Log-logistic 600    511     3 -1887.1406 3780.281
#> 9   dependent Both arms  Exponential 600    511     2 -1910.0074 3824.015
#>        BIC
#> 1 1874.541
#> 2 1884.421
#> 3 1904.804
#> 4 1907.965
#> 5 1911.605
#> 6 1926.619
#> 7 3779.105
#> 8 3793.472
#> 9 3832.809
```
