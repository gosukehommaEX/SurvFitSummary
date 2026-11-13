# Reproduction of published results

One script per published article whose tables or figures are reproduced with
the functions of the package: `reproduce_<Author><Year>.R`. Each script
writes `results/reproduce_<Author><Year>.csv` with one row per published
value (published value, value of the package, difference, verdict PASS,
EXPLAINED, FAIL or INFO). The tolerance is half a unit of the last printed
digit. A value that cannot be reproduced is EXPLAINED only with a reason
registered in `../validation/explained.csv`, and only up to the largest
difference registered with the reason (`max_difference`).

| Script | Article | Values |
|---|---|---|
| `reproduce_Jackson2016.R` | Jackson CH (2016). flexsurv: a platform for parametric survival modeling in R. *Journal of Statistical Software* 70(8), 1-33. | Breast cancer data (`flexsurv::bc`): estimates, 95% CI and standard errors of the Weibull model fs1; -2 log-likelihood, parameters and AIC of Table 3 (Weibull, generalized gamma, log-logistic, spline models sp1 and sp3); effects of sp3; 5-year survival under sp1 |
| `reproduce_RoystonParmar2002.R` | Royston P, Parmar MKB (2002). Flexible parametric proportional-hazards and proportional-odds models for censored survival data. *Statistics in Medicine* 21, 2175-2197. | Same data: AIC and group coefficients of the proportional hazards and proportional odds models with 1 to 4 degrees of freedom (Table II); knots and coefficients of the PO(2) model (Section 5.1) |

The scripts are run by `../validation/run_validation.R` (normal and
mutation runs), or alone from the root of the package after
`devtools::load_all()`, for example
`source(file.path("inst", "reproduce", "reproduce_Jackson2016.R"))`.
