# Covariance Matrices and Cholesky Factors of the Fitted Models

Reports, for each fitted model, the covariance matrix of the parameter
estimates or its Cholesky factor, on the transformed (real-line) scale
on which flexsurv estimates the parameters and/or on the natural scale.

## Usage

``` r
covariance_table(
  x,
  type = c("covariance", "cholesky"),
  scale = c("transformed", "natural")
)
```

## Arguments

- x:

  Object returned by
  [`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md).

- type:

  `"covariance"` (the default) or `"cholesky"`.

- scale:

  Character vector with `"transformed"` and/or `"natural"`. The default
  is both.

## Value

A data frame with the columns `fit_type`, `arm`, `model`, `scale`,
`row`, `col` (parameter names) and `value`, with one row per element of
the matrix of each fitted model, row by row (the elements above the
diagonal of the Cholesky factor are 0). Values are `NA` when the
covariance matrix is not available or not positive definite.

## Details

The covariance matrix on the transformed scale is the inverse of the
Hessian of minus the log-likelihood computed by flexsurv. On the natural
scale it is `J V J` (delta method), with `J` the diagonal matrix of the
derivatives of the natural parameters with respect to the transformed
ones (the estimate itself for a log-transformed parameter, 1 otherwise).
The Cholesky factor is the lower triangular matrix `L` with `L %*% t(L)`
equal to the covariance matrix; correlated normal draws for
probabilistic sensitivity analysis are `est + L z` with `z` independent
standard normal. Draws should be made on the transformed scale, where
the asymptotic normal distribution is used by flexsurv and positive
parameters stay positive after back-transformation. The parameters are
those of
[`parameter_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/parameter_table.md).

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(dat, distributions = "weibull",
                            spline_scales = NULL)
covariance_table(fits)
#>       fit_type       arm   model       scale          row          col
#> 1  independent   Control Weibull transformed        shape        shape
#> 2  independent   Control Weibull transformed        shape        scale
#> 3  independent   Control Weibull transformed        scale        shape
#> 4  independent   Control Weibull transformed        scale        scale
#> 5  independent   Control Weibull     natural        shape        shape
#> 6  independent   Control Weibull     natural        shape        scale
#> 7  independent   Control Weibull     natural        scale        shape
#> 8  independent   Control Weibull     natural        scale        scale
#> 9  independent Treatment Weibull transformed        shape        shape
#> 10 independent Treatment Weibull transformed        shape        scale
#> 11 independent Treatment Weibull transformed        scale        shape
#> 12 independent Treatment Weibull transformed        scale        scale
#> 13 independent Treatment Weibull     natural        shape        shape
#> 14 independent Treatment Weibull     natural        shape        scale
#> 15 independent Treatment Weibull     natural        scale        shape
#> 16 independent Treatment Weibull     natural        scale        scale
#> 17   dependent Both arms Weibull transformed        shape        shape
#> 18   dependent Both arms Weibull transformed        shape        scale
#> 19   dependent Both arms Weibull transformed        shape ARMTreatment
#> 20   dependent Both arms Weibull transformed        scale        shape
#> 21   dependent Both arms Weibull transformed        scale        scale
#> 22   dependent Both arms Weibull transformed        scale ARMTreatment
#> 23   dependent Both arms Weibull transformed ARMTreatment        shape
#> 24   dependent Both arms Weibull transformed ARMTreatment        scale
#> 25   dependent Both arms Weibull transformed ARMTreatment ARMTreatment
#> 26   dependent Both arms Weibull     natural        shape        shape
#> 27   dependent Both arms Weibull     natural        shape        scale
#> 28   dependent Both arms Weibull     natural        shape ARMTreatment
#> 29   dependent Both arms Weibull     natural        scale        shape
#> 30   dependent Both arms Weibull     natural        scale        scale
#> 31   dependent Both arms Weibull     natural        scale ARMTreatment
#> 32   dependent Both arms Weibull     natural ARMTreatment        shape
#> 33   dependent Both arms Weibull     natural ARMTreatment        scale
#> 34   dependent Both arms Weibull     natural ARMTreatment ARMTreatment
#>            value
#> 1   0.0020457891
#> 2   0.0005688729
#> 3   0.0005688729
#> 4   0.0022048687
#> 5   0.0036884269
#> 6   0.0100030125
#> 7   0.0100030125
#> 8   0.3781234016
#> 9   0.0024826166
#> 10  0.0004223837
#> 11  0.0004223837
#> 12  0.0025347636
#> 13  0.0042000214
#> 14  0.0114167380
#> 15  0.0114167380
#> 16  1.0946249816
#> 17  0.0011233728
#> 18  0.0003071388
#> 19 -0.0001092323
#> 20  0.0003071388
#> 21  0.0021912289
#> 22 -0.0021371197
#> 23 -0.0001092323
#> 24 -0.0021371197
#> 25  0.0044973179
#> 26  0.0019671504
#> 27  0.0053011552
#> 28 -0.0001445466
#> 29  0.0053011552
#> 30  0.3727742405
#> 31 -0.0278745540
#> 32 -0.0001445466
#> 33 -0.0278745540
#> 34  0.0044973179
covariance_table(fits, type = "cholesky", scale = "transformed")
#>       fit_type       arm   model       scale          row          col
#> 1  independent   Control Weibull transformed        shape        shape
#> 2  independent   Control Weibull transformed        shape        scale
#> 3  independent   Control Weibull transformed        scale        shape
#> 4  independent   Control Weibull transformed        scale        scale
#> 5  independent Treatment Weibull transformed        shape        shape
#> 6  independent Treatment Weibull transformed        shape        scale
#> 7  independent Treatment Weibull transformed        scale        shape
#> 8  independent Treatment Weibull transformed        scale        scale
#> 9    dependent Both arms Weibull transformed        shape        shape
#> 10   dependent Both arms Weibull transformed        shape        scale
#> 11   dependent Both arms Weibull transformed        shape ARMTreatment
#> 12   dependent Both arms Weibull transformed        scale        shape
#> 13   dependent Both arms Weibull transformed        scale        scale
#> 14   dependent Both arms Weibull transformed        scale ARMTreatment
#> 15   dependent Both arms Weibull transformed ARMTreatment        shape
#> 16   dependent Both arms Weibull transformed ARMTreatment        scale
#> 17   dependent Both arms Weibull transformed ARMTreatment ARMTreatment
#>           value
#> 1   0.045230401
#> 2   0.000000000
#> 3   0.012577224
#> 4   0.045240271
#> 5   0.049825862
#> 6   0.000000000
#> 7   0.008477198
#> 8   0.049627621
#> 9   0.033516753
#> 10  0.000000000
#> 11  0.000000000
#> 12  0.009163738
#> 13  0.045904845
#> 14  0.000000000
#> 15 -0.003259035
#> 16 -0.045904845
#> 17  0.048779523
```
