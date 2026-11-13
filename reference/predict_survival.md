# Predicted Survival, Hazard or Cumulative Hazard of the Fitted Models

Predicts the survival probability, the hazard or the cumulative hazard
of each fitted model at the given times, for each arm (the arm of an
independent fit, every arm for a dependent fit), with the Kaplan-Meier
estimate of the arm for comparison. Used for the landmark survival
probabilities (for example every 4 weeks up to the time horizon) and for
the curves of the plots.

## Usage

``` r
predict_survival(
  x,
  times,
  type = c("survival", "hazard", "cumhaz"),
  conf_int = FALSE,
  n_sim = 1000,
  seed = 1
)
```

## Arguments

- x:

  Object returned by
  [`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md).

- times:

  Non-negative times in the time unit of the data.

- type:

  `"survival"` (the default), `"hazard"` or `"cumhaz"`.

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
`model`, `time`, `est`, `lcl` and `ucl` (`NA` without `conf_int`) and
`km` (for `type = "survival"`, the Kaplan-Meier estimate of the arm, and
`NA` after its largest observed time; `NA` for the other types).

## Details

The predictions are those of `summary.flexsurvreg()`. The confidence
limits, when requested, are the quantiles of the predictions for `n_sim`
parameter vectors drawn from the asymptotic normal distribution of the
estimates on the transformed scale (flexsurv's
`normboot.flexsurvreg()`); the random number generator is seeded with
`seed` for each model and arm, so that the limits do not depend on the
other models fitted, and the state of the generator of the session is
restored afterwards. At time 0 the survival is set to 1, the cumulative
hazard to 0 and the hazard to `NA` (it is infinite or undefined for some
models), without calling flexsurv.

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(dat, distributions = c("exp", "weibull"),
                            spline_scales = NULL)
predict_survival(fits, times = seq(0, 60, by = 12))
#>       fit_type       arm       model time          est lcl ucl          km
#> 1  independent   Control Exponential    0 1.0000000000  NA  NA 1.000000000
#> 2  independent   Control Exponential   12 0.3748569551  NA  NA 0.385418910
#> 3  independent   Control Exponential   24 0.1405177368  NA  NA 0.078422012
#> 4  independent   Control Exponential   36 0.0526740510  NA  NA 0.023526604
#> 5  independent   Control Exponential   48 0.0197452344  NA  NA 0.007842201
#> 6  independent   Control Exponential   60 0.0074016384  NA  NA 0.007842201
#> 7  independent   Control     Weibull    0 1.0000000000  NA  NA 1.000000000
#> 8  independent   Control     Weibull   12 0.4109409831  NA  NA 0.385418910
#> 9  independent   Control     Weibull   24 0.1048148069  NA  NA 0.078422012
#> 10 independent   Control     Weibull   36 0.0204914161  NA  NA 0.023526604
#> 11 independent   Control     Weibull   48 0.0032770436  NA  NA 0.007842201
#> 12 independent   Control     Weibull   60 0.0004441375  NA  NA 0.007842201
#> 13 independent Treatment Exponential    0 1.0000000000  NA  NA 1.000000000
#> 14 independent Treatment Exponential   12 0.5508420044  NA  NA 0.607366234
#> 15 independent Treatment Exponential   24 0.3034269138  NA  NA 0.259972411
#> 16 independent Treatment Exponential   36 0.1671402894  NA  NA 0.131836509
#> 17 independent Treatment Exponential   48 0.0920678920  NA  NA 0.069565553
#> 18 independent Treatment Exponential   60 0.0507148622  NA  NA 0.046377035
#> 19 independent Treatment     Weibull    0 1.0000000000  NA  NA 1.000000000
#> 20 independent Treatment     Weibull   12 0.6128932006  NA  NA 0.607366234
#> 21 independent Treatment     Weibull   24 0.2993868071  NA  NA 0.259972411
#> 22 independent Treatment     Weibull   36 0.1295633797  NA  NA 0.131836509
#> 23 independent Treatment     Weibull   48 0.0512535877  NA  NA 0.069565553
#> 24 independent Treatment     Weibull   60 0.0188464022  NA  NA 0.046377035
#> 25   dependent   Control Exponential    0 1.0000000000  NA  NA 1.000000000
#> 26   dependent   Control Exponential   12 0.3748569558  NA  NA 0.385418910
#> 27   dependent   Control Exponential   24 0.1405177373  NA  NA 0.078422012
#> 28   dependent   Control Exponential   36 0.0526740513  NA  NA 0.023526604
#> 29   dependent   Control Exponential   48 0.0197452345  NA  NA 0.007842201
#> 30   dependent   Control Exponential   60 0.0074016385  NA  NA 0.007842201
#> 31   dependent Treatment Exponential    0 1.0000000000  NA  NA 1.000000000
#> 32   dependent Treatment Exponential   12 0.5508420051  NA  NA 0.607366234
#> 33   dependent Treatment Exponential   24 0.3034269146  NA  NA 0.259972411
#> 34   dependent Treatment Exponential   36 0.1671402901  NA  NA 0.131836509
#> 35   dependent Treatment Exponential   48 0.0920678925  NA  NA 0.069565553
#> 36   dependent Treatment Exponential   60 0.0507148625  NA  NA 0.046377035
#> 37   dependent   Control     Weibull    0 1.0000000000  NA  NA 1.000000000
#> 38   dependent   Control     Weibull   12 0.4083747280  NA  NA 0.385418910
#> 39   dependent   Control     Weibull   24 0.1063469668  NA  NA 0.078422012
#> 40   dependent   Control     Weibull   36 0.0216573560  NA  NA 0.023526604
#> 41   dependent   Control     Weibull   48 0.0036686433  NA  NA 0.007842201
#> 42   dependent   Control     Weibull   60 0.0005344304  NA  NA 0.007842201
#> 43   dependent Treatment     Weibull    0 1.0000000000  NA  NA 1.000000000
#> 44   dependent Treatment     Weibull   12 0.6177837648  NA  NA 0.607366234
#> 45   dependent Treatment     Weibull   24 0.2996369562  NA  NA 0.259972411
#> 46   dependent Treatment     Weibull   36 0.1273289844  NA  NA 0.131836509
#> 47   dependent Treatment     Weibull   48 0.0490058176  NA  NA 0.069565553
#> 48   dependent Treatment     Weibull   60 0.0173914545  NA  NA 0.046377035
```
