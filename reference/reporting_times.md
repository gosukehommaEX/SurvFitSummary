# Reporting Times of the Tables, Predictions and Plots

Returns the reporting times that
[`write_survival_workbook`](https://gosukehommaEX.github.io/SurvFitSummary/reference/write_survival_workbook.md)
and
[`write_survival_report`](https://gosukehommaEX.github.io/SurvFitSummary/reference/write_survival_report.md)
use, with their defaults filled in, so that tables, predictions and
plots made with the other functions of the package can use the same
times as the workbook and the report:

- the spacing of the reporting times of the Kaplan-Meier tables, by
  default the spacing of the time axis of the plots of the data (the
  smallest of the usual intervals of the time unit, for months 1, 2, 3,
  6, 12, 24, 60, 120 or 240, that gives at most 10 intervals up to the
  largest observed time);

- the landmark times of the predicted survival, by default 6 months and
  1, 2, 5 and 10 years in the time unit of the data (1 year = 365.25
  days = 12 months), labelled `"6 months"`, `"1 year"`, `"2 years"`,
  `"5 years"` and `"10 years"`;

- the time horizon of the extrapolation and of the restricted mean, by
  default the last landmark time.

## Usage

``` r
reporting_times(x, interval = NULL, landmarks = NULL, time_horizon = NULL)
```

## Arguments

- x:

  Object returned by
  [`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md).

- interval:

  `NULL` (the default) or the spacing of the reporting times, a positive
  number in the time unit of the data.

- landmarks:

  `NULL` (the default) or the landmark times, a vector of positive
  numbers in the time unit of the data. They are sorted and labelled
  with the time unit, for example `"36 months"`.

- time_horizon:

  `NULL` (the default) or the time horizon, a positive number in the
  time unit of the data.

## Value

A list with `interval` (number), `landmarks` (data frame with the
columns `time` and `label`), `time_horizon` (number) and `time_unit`
(the time unit of `x`, for example `"months"`).

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(dat, distributions = c("exp", "weibull"),
                            spline_scales = NULL)
times <- reporting_times(fits)
times
#> $interval
#> [1] 12
#> 
#> $landmarks
#>   time    label
#> 1    6 6 months
#> 2   12   1 year
#> 3   24  2 years
#> 4   60  5 years
#> 5  120 10 years
#> 
#> $time_horizon
#> [1] 120
#> 
#> $time_unit
#> [1] "months"
#> 
km_table(dat, interval = times$interval)
#>          ARM time n_risk n_event n_censor        surv          se         lcl
#> 1    Control    0    300       0        0 1.000000000 0.000000000 1.000000000
#> 2    Control   12    104     173       23 0.385418910 0.029195580 0.332241811
#> 3    Control   24     20      80        4 0.078422012 0.016689801 0.051675533
#> 4    Control   36      6      14        0 0.023526604 0.009468082 0.010690559
#> 5    Control   48      2       4        0 0.007842201 0.005519108 0.001974184
#> 6    Control   60      2       0        0 0.007842201 0.005519108 0.001974184
#> 7  Treatment    0    300       0        0 1.000000000 0.000000000 1.000000000
#> 8  Treatment   12    159     108       33 0.607366234 0.029624510 0.551992247
#> 9  Treatment   24     63      87        9 0.259972411 0.027607189 0.211123131
#> 10 Treatment   36     26      30        7 0.131836509 0.021872105 0.095239826
#> 11 Treatment   48     10      11        5 0.069565553 0.018076397 0.041803449
#> 12 Treatment   60      4       3        3 0.046377035 0.016452659 0.023138318
#>           ucl
#> 1  1.00000000
#> 2  0.44710729
#> 3  0.11901207
#> 4  0.05177476
#> 5  0.03115217
#> 6  0.03115217
#> 7  1.00000000
#> 8  0.66829515
#> 9  0.32012435
#> 10 0.18249577
#> 11 0.11576476
#> 12 0.09295531
predict_survival(fits, times = times$landmarks$time)
#>       fit_type       arm       model time          est lcl ucl          km
#> 1  independent   Control Exponential    6 6.122556e-01  NA  NA 0.729690907
#> 2  independent   Control Exponential   12 3.748570e-01  NA  NA 0.385418910
#> 3  independent   Control Exponential   24 1.405177e-01  NA  NA 0.078422012
#> 4  independent   Control Exponential   60 7.401638e-03  NA  NA 0.007842201
#> 5  independent   Control Exponential  120 5.478425e-05  NA  NA          NA
#> 6  independent   Control     Weibull    6 7.042451e-01  NA  NA 0.729690907
#> 7  independent   Control     Weibull   12 4.109410e-01  NA  NA 0.385418910
#> 8  independent   Control     Weibull   24 1.048148e-01  NA  NA 0.078422012
#> 9  independent   Control     Weibull   60 4.441375e-04  NA  NA 0.007842201
#> 10 independent   Control     Weibull  120 3.140832e-09  NA  NA          NA
#> 11 independent Treatment Exponential    6 7.421873e-01  NA  NA 0.836248773
#> 12 independent Treatment Exponential   12 5.508420e-01  NA  NA 0.607366234
#> 13 independent Treatment Exponential   24 3.034269e-01  NA  NA 0.259972411
#> 14 independent Treatment Exponential   60 5.071486e-02  NA  NA 0.046377035
#> 15 independent Treatment Exponential  120 2.571997e-03  NA  NA          NA
#> 16 independent Treatment     Weibull    6 8.197703e-01  NA  NA 0.836248773
#> 17 independent Treatment     Weibull   12 6.128932e-01  NA  NA 0.607366234
#> 18 independent Treatment     Weibull   24 2.993868e-01  NA  NA 0.259972411
#> 19 independent Treatment     Weibull   60 1.884640e-02  NA  NA 0.046377035
#> 20 independent Treatment     Weibull  120 5.637779e-05  NA  NA          NA
#> 21   dependent   Control Exponential    6 6.122556e-01  NA  NA 0.729690907
#> 22   dependent   Control Exponential   12 3.748570e-01  NA  NA 0.385418910
#> 23   dependent   Control Exponential   24 1.405177e-01  NA  NA 0.078422012
#> 24   dependent   Control Exponential   60 7.401638e-03  NA  NA 0.007842201
#> 25   dependent   Control Exponential  120 5.478425e-05  NA  NA          NA
#> 26   dependent Treatment Exponential    6 7.421873e-01  NA  NA 0.836248773
#> 27   dependent Treatment Exponential   12 5.508420e-01  NA  NA 0.607366234
#> 28   dependent Treatment Exponential   24 3.034269e-01  NA  NA 0.259972411
#> 29   dependent Treatment Exponential   60 5.071486e-02  NA  NA 0.046377035
#> 30   dependent Treatment Exponential  120 2.571997e-03  NA  NA          NA
#> 31   dependent   Control     Weibull    6 6.991509e-01  NA  NA 0.729690907
#> 32   dependent   Control     Weibull   12 4.083747e-01  NA  NA 0.385418910
#> 33   dependent   Control     Weibull   24 1.063470e-01  NA  NA 0.078422012
#> 34   dependent   Control     Weibull   60 5.344304e-04  NA  NA 0.007842201
#> 35   dependent   Control     Weibull  120 6.485911e-09  NA  NA          NA
#> 36   dependent Treatment     Weibull    6 8.249238e-01  NA  NA 0.836248773
#> 37   dependent Treatment     Weibull   12 6.177838e-01  NA  NA 0.607366234
#> 38   dependent Treatment     Weibull   24 2.996370e-01  NA  NA 0.259972411
#> 39   dependent Treatment     Weibull   60 1.739145e-02  NA  NA 0.046377035
#> 40   dependent Treatment     Weibull  120 3.950647e-05  NA  NA          NA
reporting_times(fits, interval = 6, landmarks = c(24, 36),
                time_horizon = 60)
#> $interval
#> [1] 6
#> 
#> $landmarks
#>   time     label
#> 1   24 24 months
#> 2   36 36 months
#> 
#> $time_horizon
#> [1] 60
#> 
#> $time_unit
#> [1] "months"
#> 
```
