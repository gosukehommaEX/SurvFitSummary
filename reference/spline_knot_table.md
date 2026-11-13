# Placement of the Knots of the Spline Models

Reports the internal knots of the Royston-Parmar spline models of each
fit (each arm of the independent fits and both arms together for the
dependent fit) with the Kaplan-Meier estimate at each knot, as in the
table of knot placement of a survival analysis report.

## Usage

``` r
spline_knot_table(x)
```

## Arguments

- x:

  Object returned by
  [`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md).

## Value

A data frame with the columns `fit_type`, `arm`, `n_knots`, `knot`
(index of the knot), `quantile` (probability of the quantile of the log
event times; for knots given by the user, the proportion of log event
times at or below the knot), `log_time`, `time`, `n_risk`, `surv_est`
and `surv_se`. It has no rows when no spline model was fitted.

## Details

The knots are the same for the hazard, odds and normal scales. The
number at risk and the Kaplan-Meier estimate are those of the data
fitted (both arms pooled for the dependent fit) at the knot time; a knot
at an observed time counts the patients with that time as at risk and
their events in the estimate. The standard error is that of Greenwood.

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(dat, distributions = NULL,
                            spline_scales = "hazard")
spline_knot_table(fits)
#>       fit_type       arm n_knots knot  quantile log_time      time n_risk
#> 1  independent   Control       1    1 0.5000000 2.243398  9.425300    145
#> 2  independent   Control       2    1 0.3333333 1.969362  7.166100    191
#> 3  independent   Control       2    2 0.6666667 2.542782 12.715000     97
#> 4  independent   Control       3    1 0.2500000 1.615682  5.031319    217
#> 5  independent   Control       3    2 0.5000000 2.243398  9.425300    145
#> 6  independent   Control       3    3 0.7500000 2.747653 15.605967     72
#> 7  independent Treatment       1    1 0.5000000 2.587478 13.296201    147
#> 8  independent Treatment       2    1 0.3333333 2.288967  9.864741    190
#> 9  independent Treatment       2    2 0.6666667 2.861626 17.489943    100
#> 10 independent Treatment       3    1 0.2500000 2.080707  8.010130    212
#> 11 independent Treatment       3    2 0.5000000 2.587478 13.296201    147
#> 12 independent Treatment       3    3 0.7500000 3.067678 21.491936     78
#> 13   dependent Both arms       1    1 0.5000000 2.411520 11.150900    292
#> 14   dependent Both arms       2    1 0.3333333 2.085771  8.050800    382
#> 15   dependent Both arms       2    2 0.6666667 2.727997 15.302200    198
#> 16   dependent Both arms       3    1 0.2500000 1.838845  6.289272    428
#> 17   dependent Both arms       3    2 0.5000000 2.411520 11.150900    292
#> 18   dependent Both arms       3    3 0.7500000 2.860528 17.470746    150
#>     surv_est    surv_se
#> 1  0.5205484 0.02978855
#> 2  0.6830815 0.02752297
#> 3  0.3557713 0.02876996
#> 4  0.7650842 0.02495343
#> 5  0.5205484 0.02978855
#> 6  0.2731327 0.02695286
#> 7  0.5615273 0.03019978
#> 8  0.7134789 0.02717195
#> 9  0.4037366 0.03033069
#> 10 0.7881927 0.02434523
#> 11 0.5615273 0.03019978
#> 12 0.3218706 0.02918995
#> 13 0.5421976 0.02118619
#> 14 0.6990478 0.01930597
#> 15 0.3808919 0.02090650
#> 16 0.7775871 0.01736614
#> 17 0.5421976 0.02118619
#> 18 0.2985694 0.01987936
```
