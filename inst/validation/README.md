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
| `validate_python_reference.R` | The 12 distributions and the spline models with 1 to 3 knots on the three scales, fitted to `inst/extdata/fit_example.csv`, against the independent implementation `python/expected_fits.py`: fit statistics, parameters, covariance matrices, survival and hazard, median, mean and restricted mean, knots (the tolerances of the tests, capped at half the shift of the mutation run) |
| `validate_coverage.R` | Coverage of the 95% confidence intervals of the parameters and of the predicted survival in 500 simulated samples of each of three distributions (INFO). Run alone: it takes several minutes and is not run by `run_validation.R` |
| `validate_outputs.R` | Writes the workbook and the report of `inst/extdata/fit_example.csv` (all models of `fit_survival_models()` with its defaults) to `results/outputs/`, for the checks below. Run alone (2 to 3 minutes) |

`run_validation.R` runs the scripts above except `validate_coverage.R` and
`validate_outputs.R`, and those of `../reproduce/`, twice:
normally, where no comparison may FAIL and every entry of `explained.csv`
must be used, and with the values of the package shifted by
0.01 * max(1, |value|) (environment variable
`SURVFITSUMMARY_MUTATION = "1"`, results in `results/mutation/`), where no
comparison may PASS or be EXPLAINED. Every tolerance must therefore be at
most 0.005 * max(1, |reference|), half the shift, and `verdict.R` stops a
script whose tolerance is wider. It writes the summary
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

`python/check_workbook.py` and `python/check_report.py` check the workbook
and the report written by `validate_outputs.R` (numbers against the data
and `python/expected/`, numbering of the tables and figures; openpyxl and
python-docx). From the root of the package:

```
python3 -I inst/validation/python/check_workbook.py inst/validation/results/outputs/fit_example_results.xlsx inst/validation/python/expected inst/extdata/fit_example.csv --csv inst/validation/results/outputs/check_workbook.csv
python3 -I inst/validation/python/check_workbook.py inst/validation/results/outputs/fit_example_results.xlsx inst/validation/python/expected inst/extdata/fit_example.csv --mutate --csv inst/validation/results/outputs/mutation/check_workbook.csv
```

and likewise `check_report.py` with `fit_example_report.docx` and
`check_report.csv`. With `--mutate`, every number read from the file (except
in the columns that identify a row) is shifted as in the mutation run, and
every check must FAIL; the checks of the numbering are left out. The
workbook and the report are not kept in git or in the package
(`.gitignore`, `.Rbuildignore`; the report is about 16 MB); the verdicts
are.
