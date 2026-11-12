# Validation by independent routes

One script per comparison of independent routes (closed form against
simulation, the package against another package, the package against an
independent implementation): `validate_<content>.R`. Each script writes
`results/validate_<content>.csv` with the verdict of each comparison (PASS,
EXPLAINED, FAIL or INFO). The reasons of the EXPLAINED verdicts are
registered in `explained.csv`; the summary of the verdicts is written by the
scripts, never by hand.

`python/` holds the independent implementations (numpy and scipy, no
survival analysis library) that give the expected values of the tests:
`expected_km_cox.py` (Kaplan-Meier, Cox and Grambsch-Therneau tests; prints
the values written in the tests), `make_fit_example.py` (simulates
`inst/extdata/fit_example.csv`) and `expected_fits.py` (maximum likelihood
fits of the parametric and spline models; writes `python/expected/*.csv`,
read by the tests). Run them with `python3 -I <script>`.
