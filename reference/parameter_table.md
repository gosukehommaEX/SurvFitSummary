# Parameter Estimates of the Fitted Models

Reports the maximum likelihood estimates of the parameters of each
fitted model with their standard errors and Wald confidence intervals,
on the natural scale of the distribution and/or on the transformed
(real-line) scale on which flexsurv estimates them.

## Usage

``` r
parameter_table(x, scale = c("natural", "transformed"))
```

## Arguments

- x:

  Object returned by
  [`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md).

- scale:

  Character vector with `"natural"` and/or `"transformed"`. The default
  is both.

## Value

A data frame with the columns `fit_type`, `arm`, `model`, `parameter`
(name used by flexsurv, for example `"shape"`, `"gamma1"` or
`"ARMTreatment"`), `scale`, `est`, `lcl`, `ucl` and `se`, with one row
per parameter and scale of each fitted model (models that failed are
left out).

## Details

On the transformed scale, the positive parameters (for example the shape
and scale of the Weibull, see
[`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md))
are on the log scale and the others, the arm effects and all spline
parameters are untransformed; the confidence interval is
`est +/- z * se` with `z` the normal quantile of the confidence level of
the fit. On the natural scale, the estimate and the confidence limits of
a log-transformed parameter are the exponentials of those on the
transformed scale (as in `flexsurvreg()$res`) and the standard error is
obtained by the delta method (`est * se` of the log). Arm effects in
dependent fits are effects on the transformed location parameter (for
example the log of the scale of the Weibull) on both scales. The shape
of the Rayleigh model, fixed at 2, is not reported.

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(dat, distributions = c("exp", "weibull"),
                            spline_scales = NULL)
parameter_table(fits)
#>       fit_type       arm       model    parameter       scale         est
#> 1  independent   Control Exponential         rate     natural  0.08176756
#> 2  independent   Control Exponential         rate transformed -2.50387463
#> 3  independent   Control     Weibull        shape     natural  1.34273452
#> 4  independent   Control     Weibull        scale     natural 13.09560028
#> 5  independent   Control     Weibull        shape transformed  0.29470822
#> 6  independent   Control     Weibull        scale transformed  2.57227632
#> 7  independent Treatment Exponential         rate     natural  0.04969227
#> 8  independent Treatment Exponential         rate transformed -3.00190587
#> 9  independent Treatment     Weibull        shape     natural  1.30068139
#> 10 independent Treatment     Weibull        scale     natural 20.78088053
#> 11 independent Treatment     Weibull        shape transformed  0.26288827
#> 12 independent Treatment     Weibull        scale transformed  3.03403336
#> 13   dependent Both arms Exponential         rate     natural  0.08176756
#> 14   dependent Both arms Exponential ARMTreatment     natural -0.49803124
#> 15   dependent Both arms Exponential         rate transformed -2.50387463
#> 16   dependent Both arms Exponential ARMTreatment transformed -0.49803124
#> 17   dependent Both arms     Weibull        shape     natural  1.32329549
#> 18   dependent Both arms     Weibull        scale     natural 13.04304754
#> 19   dependent Both arms     Weibull ARMTreatment     natural  0.46876283
#> 20   dependent Both arms     Weibull        shape transformed  0.28012521
#> 21   dependent Both arms     Weibull        scale transformed  2.56825524
#> 22   dependent Both arms     Weibull ARMTreatment transformed  0.46876283
#>            lcl         ucl          se
#> 1   0.07258957  0.09210599 0.004967026
#> 2  -2.62293396 -2.38481530 0.060745674
#> 3   1.22882480  1.46720346 0.060732421
#> 4  11.94418071 14.35801674 0.614917394
#> 5   0.20605827  0.38335818 0.045230401
#> 6   2.48024419  2.66430844 0.046956030
#> 7   0.04378688  0.05639410 0.003207622
#> 8  -3.12842100 -2.87539074 0.064549722
#> 9   1.17966603  1.43411103 0.064807572
#> 10 18.82820775 22.93606494 1.046243271
#> 11  0.16523138  0.36054517 0.049825862
#> 12  2.93535616  3.13271056 0.050346436
#> 13  0.07258957  0.09210599 0.004967026
#> 14 -0.67175862 -0.32430385 0.088638048
#> 15 -2.62293397 -2.38481530 0.060745674
#> 16 -0.67175862 -0.32430385 0.088638048
#> 17  1.23915981  1.41314376 0.044352569
#> 18 11.89964083 14.29632133 0.610552406
#> 19  0.33732364  0.60020202 0.067062045
#> 20  0.21443358  0.34581684 0.033516753
#> 21  2.47650822  2.66000225 0.046810564
#> 22  0.33732364  0.60020202 0.067062045
```
