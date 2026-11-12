# Validation by independent routes

One script per comparison of independent routes (closed form against
simulation, the package against another package, the package against an
independent implementation): `validate_<content>.R`. Each script writes
`results/validate_<content>.csv` with the verdict of each comparison (PASS,
EXPLAINED, FAIL or INFO). The reasons of the EXPLAINED verdicts are
registered in `explained.csv`; the summary of the verdicts is written by the
scripts, never by hand.
