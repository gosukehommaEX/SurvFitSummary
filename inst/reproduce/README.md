# Reproduction of published results

One script per published article whose tables or figures are reproduced with
the functions of the package: `reproduce_<Author><Year>.R`. Each script
writes `results/reproduce_<Author><Year>.csv` with one row per published
value (published value, value of the package, difference, verdict PASS,
EXPLAINED, FAIL or INFO). A value that cannot be reproduced is EXPLAINED only
with a reason registered in `../validation/explained.csv`.
