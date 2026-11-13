# Validation of SurvFitSummary

``` r

library(SurvFitSummary)
read_result <- function(...) {
  utils::read.csv(system.file(..., package = "SurvFitSummary",
                              mustWork = TRUE),
                  stringsAsFactors = FALSE)
}
f1 <- function(x) sprintf("%.1f", x)
f3 <- function(x) sprintf("%.3f", x)
runs <- read_result("validation", "results", "validation_summary.csv")
normal <- runs[runs$run == "normal", ]
mutation <- runs[runs$run == "mutation", ]
```

The results of SurvFitSummary are checked in three ways:

1.  The unit tests (`tests/testthat`) compare every fitted model, table
    and output file with reference values computed independently,
    without any survival analysis library, by the Python scripts in
    `system.file("validation", "python", package = "SurvFitSummary")`.
2.  The validation scripts (`inst/validation`) compare the package with
    closed forms, with other R packages that implement the same models
    independently, and with the Python reference, and estimate the
    coverage of the confidence intervals by simulation.
3.  The reproduction scripts (`inst/reproduce`) reproduce published
    results: those of Jackson (2016) and of Royston and Parmar (2002)
    for the breast cancer data
    [`flexsurv::bc`](http://chjackson.github.io/flexsurv-dev/reference/bc.md).

The validation and reproduction scripts write the verdict of every
comparison to the CSV files that are installed with the package and read
by this vignette. The latest run gave 5149 PASS, 7 EXPLAINED, 0 FAIL and
3 INFO verdicts.

## Verdicts

Each comparison is PASS if the difference from the reference is within
its tolerance, EXPLAINED if it is not but the difference is registered
with its reason and an upper bound in `explained.csv`, FAIL otherwise,
and INFO for a value reported without a reference. The tolerances of the
published values are half a unit of their last printed digit.

``` r

knitr::kable(normal[, c("script", "PASS", "EXPLAINED", "FAIL", "INFO")],
             row.names = FALSE)
```

| script                      | PASS | EXPLAINED | FAIL | INFO |
|:----------------------------|-----:|----------:|-----:|-----:|
| validate_survreg            |  150 |         0 |    0 |    0 |
| validate_closed_forms       |  186 |         0 |    0 |    2 |
| validate_splines_rstpm2     |  312 |         0 |    0 |    0 |
| validate_python_reference   | 4429 |         0 |    0 |    0 |
| reproduce_Jackson2016       |   35 |         6 |    0 |    0 |
| reproduce_RoystonParmar2002 |   37 |         1 |    0 |    1 |

## Independent routes

| Script | Comparison |
|----|----|
| `validate_survreg.R` | Exponential, Weibull, lognormal and log-logistic models against [`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html) |
| `validate_closed_forms.R` | Exponential maximum likelihood estimates in closed form; survival, hazard, median, mean and restricted mean of six distributions against their formulas |
| `validate_splines_rstpm2.R` | Spline models on the three scales against [`rstpm2::stpm2()`](https://rdrr.io/pkg/rstpm2/man/gsm.html) with the same knots |
| `validate_python_reference.R` | The 12 distributions and the spline models against the independent Python implementation |

The largest difference of the comparisons that pass, as a fraction of
their tolerance, shows how much margin each comparison has.

``` r

scripts <- c("validate_survreg", "validate_closed_forms",
             "validate_splines_rstpm2", "validate_python_reference")
margins <- do.call(rbind, lapply(scripts, function(s) {
  v <- read_result("validation", "results", paste0(s, ".csv"))
  p <- v[which(v$verdict == "PASS" & v$tolerance > 0), ]
  data.frame(script = s, comparisons = nrow(v),
             largest_fraction_of_tolerance = max(p$difference / p$tolerance))
}))
knitr::kable(margins, digits = 3)
```

| script                    | comparisons | largest_fraction_of_tolerance |
|:--------------------------|------------:|------------------------------:|
| validate_survreg          |         150 |                         0.001 |
| validate_closed_forms     |         188 |                         0.006 |
| validate_splines_rstpm2   |         312 |                         0.090 |
| validate_python_reference |        4429 |                         0.181 |

## Published results

``` r

repro <- rbind(
  read_result("reproduce", "results", "reproduce_Jackson2016.csv"),
  read_result("reproduce", "results", "reproduce_RoystonParmar2002.csv")
)
compared <- repro[repro$verdict != "INFO", ]
explained <- read_result("validation", "explained.csv")
```

Of the 79 published values, 72 are reproduced to their printed digits
and 7 differ for the reasons registered in `explained.csv`: the
published fits of the spline models of Jackson (2016) differ from the
maximum likelihood estimates, on which the package and `rstpm2` agree,
in the third decimal of the effects and the sixth decimal of the
survival probabilities (the article computes the 5-year survival of the
odds-scale model on the hazard scale, which the reference value takes
into account), and one estimate of Royston and Parmar (2002) lies just
beyond the rounding boundary of its printed value.

``` r

ex <- compared[compared$verdict == "EXPLAINED",
               c("script", "item", "reference", "value", "difference",
                 "max_difference")]
knitr::kable(ex, digits = 6, row.names = FALSE)
```

| script | item | reference | value | difference | max_difference |
|:---|:---|---:|---:|---:|---:|
| reproduce_Jackson2016 | sp3: ARMMedium est | 0.833452 | 0.834503 | 0.001051 | 2e-03 |
| reproduce_Jackson2016 | sp3: ARMMedium se | 0.171204 | 0.171276 | 0.000072 | 2e-04 |
| reproduce_Jackson2016 | sp3: ARMPoor est | 1.611179 | 1.612099 | 0.000920 | 2e-03 |
| reproduce_Jackson2016 | sp3: ARMPoor se | 0.164093 | 0.164175 | 0.000082 | 2e-04 |
| reproduce_Jackson2016 | sp1: 5-year survival of the Good group (odds scale, from the printed value) | 0.729127 | 0.729131 | 0.000004 | 2e-05 |
| reproduce_Jackson2016 | sp1: 5-year survival as computed in the article (hazard scale) | 0.689697 | 0.689702 | 0.000006 | 2e-05 |
| reproduce_RoystonParmar2002 | Table II, PO(3) beta2 (Poor) | 2.162000 | 2.162522 | 0.000522 | 1e-03 |

``` r

knitr::kable(explained[, c("item", "reason")], row.names = FALSE)
```

| item | reason |
|:---|:---|
| sp3: ARMMedium est | The published fit (flexsurv, 2016) differs from the maximum likelihood estimate in the third decimal of the effects and the fifth decimal of their standard errors; the package (optim reltol 1e-10) and rstpm2 agree with each other. |
| sp3: ARMMedium se | The published fit (flexsurv, 2016) differs from the maximum likelihood estimate in the third decimal of the effects and the fifth decimal of their standard errors; the package (optim reltol 1e-10) and rstpm2 agree with each other. |
| sp3: ARMPoor est | The published fit (flexsurv, 2016) differs from the maximum likelihood estimate in the third decimal of the effects and the fifth decimal of their standard errors; the package (optim reltol 1e-10) and rstpm2 agree with each other. |
| sp3: ARMPoor se | The published fit (flexsurv, 2016) differs from the maximum likelihood estimate in the third decimal of the effects and the fifth decimal of their standard errors; the package (optim reltol 1e-10) and rstpm2 agree with each other. |
| sp1: 5-year survival of the Good group (odds scale, from the printed value) | The article computes 1 - psurvspline(5, gamma, knots), whose default scale is “hazard”, for the model sp1 fitted on the odds scale: the printed 0.6896969 is exp(-exp(eta)), from which the survival of sp1 is 1 / (1 - log(0.6896969)) = 0.729127. The remaining difference comes from the published fit, as for sp3. |
| sp1: 5-year survival as computed in the article (hazard scale) | Same computation as the article; the remaining difference comes from the published fit, as for sp3. |
| Table II, PO(3) beta2 (Poor) | The estimate 2.16252 is 2e-5 beyond the rounding boundary 2.1625 of the printed 2.162: the fit of the article (Stata) and that of the package differ at least in the fifth decimal (possible causes: the optimizer, and the knots at the 33.3% and 66.7% centiles of the log event times, whose definition may differ between Stata and quantile(type = 7) of R). |

## Mutation run

Every script is run a second time with the values of the package shifted
by 0.01 \* max(1, \|value\|) before the comparison. Every comparison
must then FAIL; a comparison that still passed would not detect an error
of that size. For this, every tolerance is at most half the shift, which
the scripts check before they compare. In the latest mutation run there
were 0 PASS and 0 EXPLAINED verdicts among the 5156 comparisons.

``` r

knitr::kable(mutation[, c("script", "PASS", "EXPLAINED", "FAIL", "INFO")],
             row.names = FALSE)
```

| script                      | PASS | EXPLAINED | FAIL | INFO |
|:----------------------------|-----:|----------:|-----:|-----:|
| validate_survreg            |    0 |         0 |  150 |    0 |
| validate_closed_forms       |    0 |         0 |  186 |    2 |
| validate_splines_rstpm2     |    0 |         0 |  312 |    0 |
| validate_python_reference   |    0 |         0 | 4429 |    0 |
| reproduce_Jackson2016       |    0 |         0 |   41 |    0 |
| reproduce_RoystonParmar2002 |    0 |         0 |   38 |    1 |

## Coverage of the confidence intervals

`validate_coverage.R` simulates 500 samples of 300 patients from
Weibull, lognormal and log-logistic distributions, fits each with its
own distribution (and the Weibull samples also with the spline model on
the hazard scale with 1 knot), and records whether the 95% confidence
intervals of the parameters and of the survival at 6, 12 and 24 months
contain the true values.

``` r

cov <- read_result("validation", "results", "validate_coverage.csv")
coverage <- cov[grepl("coverage of the 95% CI", cov$item), ]
mc_se <- cov$value[grepl("Monte Carlo", cov$item)]
knitr::kable(data.frame(interval = coverage$item,
                        coverage = coverage$value),
             digits = 3)
```

| interval | coverage |
|:---|---:|
| Weibull samples, Weibull model: coverage of the 95% CI of shape | 0.960 |
| Weibull samples, Weibull model: coverage of the 95% CI of scale | 0.950 |
| Weibull samples, Weibull model: coverage of the 95% CI of S(6) | 0.956 |
| Weibull samples, Weibull model: coverage of the 95% CI of S(12) | 0.948 |
| Weibull samples, Weibull model: coverage of the 95% CI of S(24) | 0.954 |
| Weibull samples, Hazards, 1 knot model: coverage of the 95% CI of gamma0 | 0.964 |
| Weibull samples, Hazards, 1 knot model: coverage of the 95% CI of gamma1 | 0.936 |
| Weibull samples, Hazards, 1 knot model: coverage of the 95% CI of gamma2 | 0.946 |
| Weibull samples, Hazards, 1 knot model: coverage of the 95% CI of S(6) | 0.964 |
| Weibull samples, Hazards, 1 knot model: coverage of the 95% CI of S(12) | 0.944 |
| Weibull samples, Hazards, 1 knot model: coverage of the 95% CI of S(24) | 0.948 |
| Lognormal samples, Lognormal model: coverage of the 95% CI of meanlog | 0.940 |
| Lognormal samples, Lognormal model: coverage of the 95% CI of sdlog | 0.950 |
| Lognormal samples, Lognormal model: coverage of the 95% CI of S(6) | 0.946 |
| Lognormal samples, Lognormal model: coverage of the 95% CI of S(12) | 0.936 |
| Lognormal samples, Lognormal model: coverage of the 95% CI of S(24) | 0.950 |
| Log-logistic samples, Log-logistic model: coverage of the 95% CI of shape | 0.940 |
| Log-logistic samples, Log-logistic model: coverage of the 95% CI of scale | 0.944 |
| Log-logistic samples, Log-logistic model: coverage of the 95% CI of S(6) | 0.944 |
| Log-logistic samples, Log-logistic model: coverage of the 95% CI of S(12) | 0.948 |
| Log-logistic samples, Log-logistic model: coverage of the 95% CI of S(24) | 0.952 |

The coverage ranges from 0.936 to 0.964, within 1.4 Monte Carlo standard
errors (0.010) of 0.95.

## Excel workbook and Word report

`validate_outputs.R` writes the workbook and the report of the example
data with all the default models, and the Python scripts
`check_workbook.py` and `check_report.py` check their numbers against
the data and the reference values, and the numbering of the tables and
figures of the report.

19 of the 19 checks of the workbook and 10 of the 10 checks of the
report pass. With every number of the files shifted as in the mutation
run, 19 of the 19 checks of the workbook and 6 of the 6 checks of the
numbers of the report fail.

## Running the validation

The scripts are run from the root of the source package after
`devtools::load_all()`; `validate_splines_rstpm2.R` needs the package
rstpm2, which SurvFitSummary itself does not use.

``` r

source(file.path("inst", "validation", "run_validation.R"))
source(file.path("inst", "validation", "validate_coverage.R"))
source(file.path("inst", "validation", "validate_outputs.R"))
```

## References

Jackson C (2016). flexsurv: A platform for parametric survival modeling
in R. Journal of Statistical Software 70(8):1-33.

Royston P, Parmar MKB (2002). Flexible parametric proportional-hazards
and proportional-odds models for censored survival data, with
application to prognostic modelling and estimation of treatment effects.
Statistics in Medicine 21(15):2175-2197.
