# Fit Standard Parametric and Spline Survival Models

Fits parametric survival models with
[`flexsurv::flexsurvreg()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
and Royston-Parmar spline models with
[`flexsurv::flexsurvspline()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)
to each treatment arm separately (independent fits) and to all arms
together with the arm as a covariate (dependent fits, labelled
`"Both arms"`). A model that fails is kept in the list of models with
its error message and the other models are still fitted; warnings of
flexsurv are recorded in the same way instead of being printed. The
returned object is used by
[`parameter_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/parameter_table.md),
[`covariance_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/covariance_table.md),
[`gof_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/gof_table.md),
[`predict_survival`](https://gosukehommaEX.github.io/SurvFitSummary/reference/predict_survival.md),
[`survival_summary`](https://gosukehommaEX.github.io/SurvFitSummary/reference/survival_summary.md)
and
[`spline_knot_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/spline_knot_table.md).

## Usage

``` r
fit_survival_models(
  data,
  distributions = c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma"),
  spline_scales = c("hazard", "odds", "normal"),
  spline_knots = 1:3,
  knots = NULL,
  fit_type = c("independent", "dependent"),
  time_unit = c("months", "days", "weeks", "years"),
  conf_level = 0.95,
  control = list(reltol = 1e-10, maxit = 1000)
)
```

## Arguments

- data:

  Data frame returned by
  [`processing_dataset`](https://gosukehommaEX.github.io/SurvFitSummary/reference/processing_dataset.md).

- distributions:

  Character vector of distribution codes (see Details), or `NULL` for no
  parametric model. The default is the seven standard models
  `c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma")`.

- spline_scales:

  Character vector with the scales of the spline models, a subset of
  `c("hazard", "odds", "normal")` (the default), or `NULL` for no spline
  model.

- spline_knots:

  Numbers of internal knots of the spline models, whole numbers from 1
  to 3. The default is `1:3`. Ignored when `knots` is given.

- knots:

  `NULL` (the default) for knots at quantiles of the log event times, or
  a list of increasing vectors of internal knots on the time scale, for
  example `list(12, c(6, 18))`; each vector gives one spline model per
  scale, and the vectors must have different lengths. The knots must lie
  strictly between the smallest and largest event times of each fit.

- fit_type:

  Character vector with `"independent"` (each arm separately) and/or
  `"dependent"` (all arms with the arm as a covariate). The default is
  both. With a single arm only independent models are fitted, with a
  warning if dependent models were requested.

- time_unit:

  Time unit of `SURVTIME`: `"months"` (the default), `"days"`, `"weeks"`
  or `"years"`. It is stored in the result for the labels and time
  conversions of the outputs.

- conf_level:

  Confidence level, a single number in (0, 1). The default is 0.95.

- control:

  List of settings passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html). The default is
  `list(reltol = 1e-10, maxit = 1000)`; `ndeps` is added for the spline
  models on the normal scale (see Details).

## Value

An object of class `survfitsummary`, a list with

- models:

  Data frame with one row per model: `model_id`, `fit_type`
  (`"independent"` or `"dependent"`), `arm` (the arm, or `"Both arms"`
  for dependent fits), `family` (`"parametric"` or `"spline"`),
  `distribution` (code, `NA` for splines), `scale` and `n_knots` (`NA`
  for parametric models), `model` (print-ready label such as
  `"Log-logistic"` or `"Hazards, 2 knots"`), `n` and `events` (of the
  data fitted), `status` (`"ok"`, `"warning"` or `"error"`) and
  `message` (warnings or error, `NA` if none).

- fits:

  List of flexsurvreg objects in the order of `models`, `NULL` for
  models that failed.

- knots:

  Data frame of the internal knots of the spline models: `fit_type`,
  `arm`, `n_knots`, `knot` (index), `quantile` and `log_time`.

- boundary_knots:

  Data frame with `fit_type`, `arm`, `lower` and `upper` (log time).

- data:

  The columns `ARM`, `SURVTIME` and `EVENT` of `data`.

- arms:

  Treatment arms, the control arm first.

- time_unit, conf_level:

  The arguments.

## Details

The distributions, with their parameters in the order and on the natural
scale of flexsurv, are:

- `"exp"` (Exponential): rate.

- `"weibull"` (Weibull, accelerated failure time form): shape, scale;
  `S(t) = exp(-(t / scale)^shape)`.

- `"weibullPH"` (Weibull, proportional hazards form): shape, scale;
  `S(t) = exp(-scale * t^shape)`.

- `"rayleigh"` (Rayleigh): scale; the Weibull with the shape fixed at 2,
  so that `S(t) = exp(-(t / scale)^2)`. The fixed shape is not counted
  as a parameter.

- `"lnorm"` (Lognormal): meanlog, sdlog.

- `"llogis"` (Log-logistic): shape, scale;
  `S(t) = 1 / (1 + (t / scale)^shape)`.

- `"gompertz"` (Gompertz): shape, rate; hazard `rate * exp(shape * t)`.

- `"gamma"` (Gamma): shape, rate.

- `"gengamma"` (Generalized Gamma, Prentice 1974): mu, sigma, Q.

- `"gengamma.orig"` (Generalized Gamma, Stacy 1962): shape, scale, k.
  Equivalent to `"gengamma"` with Q \> 0 only.

- `"genf"` (Generalized F, Prentice 1975): mu, sigma, Q, P.

- `"genf.orig"` (Generalized F, original form): mu, sigma, s1, s2.

The default is the seven standard parametric models of NICE DSU
Technical Support Document 14. In dependent fits the arm is a covariate
of the location parameter on its transformed scale: the rate of the
exponential, Gompertz and gamma, the scale of the Weibull, Weibull (PH),
Rayleigh, log-logistic and original generalized gamma, meanlog of the
lognormal, and mu of the generalized gamma and generalized F.

The spline models (Royston and Parmar 2002) model the log cumulative
hazard (`"hazard"`), the log cumulative odds (`"odds"`) or the probit of
the survival (`"normal"`) as a natural cubic spline of log time with
parameters `gamma0`, `gamma1`, ...; in dependent fits the arm is a
covariate of `gamma0`. The boundary knots are the smallest and largest
log event times of the data fitted (both arms together for dependent
fits). By default the internal knots are at equally spaced quantiles of
the log event times: the median for one knot, the 1/3 and 2/3 quantiles
for two knots and the quartiles for three knots.

The likelihood is maximised by
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) (BFGS) with the
settings in `control`; the default tolerance is smaller than that of
[`optim()`](https://rdrr.io/r/stats/optim.html) so that the estimates
are accurate to well below their standard errors. flexsurv has no
analytic gradient for the spline models on the normal scale, for which
[`optim()`](https://rdrr.io/r/stats/optim.html) uses central
differences: their step is `ndeps = 1e-5` (instead of the default `1e-3`
of [`optim()`](https://rdrr.io/r/stats/optim.html), which stops short of
the maximum) unless `control` sets `ndeps`.

## References

Jackson C (2016). flexsurv: A platform for parametric survival modeling
in R. Journal of Statistical Software 70(8):1-33.
[doi:10.18637/jss.v070.i08](https://doi.org/10.18637/jss.v070.i08)

Latimer NR (2011, updated 2013). NICE DSU Technical Support Document 14:
Survival analysis for economic evaluations alongside clinical trials -
extrapolation with patient-level data.

Royston P, Parmar MKB (2002). Flexible parametric proportional-hazards
and proportional-odds models for censored survival data, with
application to prognostic modelling and estimation of treatment effects.
Statistics in Medicine 21(15):2175-2197.
[doi:10.1002/sim.1203](https://doi.org/10.1002/sim.1203)

## Examples

``` r
dat <- processing_dataset(
  system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(dat, distributions = c("exp", "weibull"),
                            spline_scales = "hazard", spline_knots = 1)
fits
#> Survival models fitted with SurvFitSummary (time unit: months)
#> Fits: Control (independent), Treatment (independent), Both arms (dependent)
#> 9 models: 9 without problems, 0 with warnings, 0 failed
fits$models[, c("fit_type", "arm", "model", "status")]
#>      fit_type       arm           model status
#> 1 independent   Control     Exponential     ok
#> 2 independent   Control         Weibull     ok
#> 3 independent   Control Hazards, 1 knot     ok
#> 4 independent Treatment     Exponential     ok
#> 5 independent Treatment         Weibull     ok
#> 6 independent Treatment Hazards, 1 knot     ok
#> 7   dependent Both arms     Exponential     ok
#> 8   dependent Both arms         Weibull     ok
#> 9   dependent Both arms Hazards, 1 knot     ok
```
