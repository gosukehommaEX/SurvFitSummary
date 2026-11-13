# Validation by independent routes

One script per comparison of independent routes (closed form against
simulation, the package against another package, the package against an
independent implementation): `validate_<content>.R`. Each script writes
`results/validate_<content>.csv` with the verdict of each comparison (PASS,
EXPLAINED, FAIL or INFO), using the functions of `verdict.R`. The reasons
of the EXPLAINED verdicts are registered in `explained.csv` (columns
`script`, `item`, `max_difference`, `reason` and `evidence`): an item that
misses its tolerance is EXPLAINED only if its difference is at most
`max_difference`, the largest difference the reason accounts for, which
must be below 0.005 so that the mutation run below cannot be EXPLAINED.
The verdicts are written by the scripts, never by hand.

| Script | Comparison |
|---|---|
| `validate_survreg.R` | Exponential, Weibull, lognormal and log-logistic models against `survival::survreg()` (log-likelihood, parameters and standard errors; `survival::lung` and `flexsurv::bc`, each arm and the arms together) |
| `validate_closed_forms.R` | Exponential maximum likelihood estimates in closed form; survival, hazard, median, mean and restricted mean of six distributions against their formulas |
| `validate_splines_rstpm2.R` | Spline models on the hazard, odds and normal scales against `rstpm2::stpm2()` with the same knots (log-likelihood, predicted survival, effect of the arm). Needs the package rstpm2, which the package itself does not use |

`run_validation.R` runs these scripts and those of `../reproduce/` twice:
normally, where no comparison may FAIL and every entry of `explained.csv`
must be used, and with the values of the package shifted by
0.01 * max(1, |value|) (environment variable
`SURVFITSUMMARY_MUTATION = "1"`, results in `results/mutation/`), where no
comparison may PASS or be EXPLAINED. It writes the summary
`results/validation_summary.csv`. Run it from the root of the package
after `devtools::load_all()`:

```r
source(file.path("inst", "validation", "run_validation.R"))
```

`python/` holds the independent implementations (numpy and scipy, no
survival analysis library) that give the expected values of the tests:
`expected_km_cox.py` (Kaplan-Meier, Cox and Grambsch-Therneau tests; prints
the values written in the tests), `make_fit_example.py` (simulates
`inst/extdata/fit_example.csv`) and `expected_fits.py` (maximum likelihood
fits of the parametric and spline models; writes `python/expected/*.csv`,
read by the tests). Run them with `python3 -I <script>`.

`python/check_workbook.py` and `python/check_report.py` check a workbook
and a report written for `inst/extdata/fit_example.csv` with all models
(numbers against the data and `python/expected/`, numbering of tables
and figures; openpyxl and python-docx):
`python3 -I python/check_workbook.py <file.xlsx> python/expected <fit_example.csv>`.
