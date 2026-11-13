# Summary of the Survival Models Fitted by fit_survival_models()

The [`summary()`](https://rdrr.io/r/base/summary.html) method of the
object returned by
[`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md):
the fit statistics of
[`gof_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/gof_table.md)
for all the models, with their status, including the models that failed
(whose statistics are `NA`). The warnings and errors themselves are
printed by
[`print.survfitsummary`](https://gosukehommaEX.github.io/SurvFitSummary/reference/print.survfitsummary.md).

## Usage

``` r
# S3 method for class 'survfitsummary'
summary(object, ...)
```

## Arguments

- object:

  Object returned by
  [`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md).

- ...:

  Not used.

## Value

A data frame with one row per model and the columns `fit_type`, `arm`,
`model`, `status` (`"ok"`, `"warning"` or `"error"`), `n`, `events`,
`npars`, `loglik`, `AIC` and `BIC`, sorted by AIC within each fit as in
[`gof_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/gof_table.md),
with the models that failed last.

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
summary(fits)
#>      fit_type       arm        model status   n events npars     loglik
#> 1 independent   Control      Weibull     ok 300    271     2  -931.5666
#> 2 independent   Control Log-logistic     ok 300    271     2  -936.5067
#> 3 independent   Control  Exponential     ok 300    271     1  -949.5500
#> 4 independent Treatment      Weibull     ok 300    240     2  -948.2785
#> 5 independent Treatment Log-logistic     ok 300    240     2  -950.0985
#> 6 independent Treatment  Exponential     ok 300    240     1  -960.4574
#> 7   dependent Both arms      Weibull     ok 600    511     3 -1879.9570
#> 8   dependent Both arms Log-logistic     ok 600    511     3 -1887.1406
#> 9   dependent Both arms  Exponential     ok 600    511     2 -1910.0074
#>        AIC      BIC
#> 1 1867.133 1874.541
#> 2 1877.013 1884.421
#> 3 1901.100 1904.804
#> 4 1900.557 1907.965
#> 5 1904.197 1911.605
#> 6 1922.915 1926.619
#> 7 3765.914 3779.105
#> 8 3780.281 3793.472
#> 9 3824.015 3832.809
```
