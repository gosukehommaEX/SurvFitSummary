# Smoothed Empirical Hazard by Treatment Arm

Estimates the hazard function of each treatment arm without a parametric
assumption with
[`bshazard::bshazard()`](https://rdrr.io/pkg/bshazard/man/bshazard.html),
which smooths the hazard with B-splines whose smoothness is estimated
from the data through a Poisson generalised linear mixed model (Rebora
et al., 2014).

## Usage

``` r
estimate_empirical_hazard(
  data,
  conf_level = 0.95,
  nbin = NULL,
  nk = 31,
  degree = 1
)
```

## Arguments

- data:

  Data frame returned by
  [`processing_dataset`](https://gosukehommaEX.github.io/SurvFitSummary/reference/processing_dataset.md).

- conf_level:

  Confidence level of the pointwise intervals, a single number in (0,
  1). The default is 0.95.

- nbin:

  Number of time intervals in which the data are split, passed to
  [`bshazard::bshazard()`](https://rdrr.io/pkg/bshazard/man/bshazard.html).
  The default `NULL` uses the distinct observed times.

- nk:

  Number of knots of the B-splines. The default is 31, as in
  [`bshazard::bshazard()`](https://rdrr.io/pkg/bshazard/man/bshazard.html).

- degree:

  Degree of the B-splines. The default is 1, as in
  [`bshazard::bshazard()`](https://rdrr.io/pkg/bshazard/man/bshazard.html).

## Value

A data frame with the columns `ARM`, `time`, `hazard`, `lcl` and `ucl`.
Rows with a non-finite hazard are removed.

## References

Rebora P, Salim A, Reilly M (2014). bshazard: A flexible tool for
nonparametric smoothing of the hazard function. The R Journal
6(2):114-122.
[doi:10.32614/RJ-2014-028](https://doi.org/10.32614/RJ-2014-028)

## Examples

``` r
raw <- generate_dummy_survival_data(
  arm = c("Control", "Treatment"), n = c(100, 100),
  hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 123
)
dat <- processing_dataset(raw, column_arm = "ARM",
                          column_survtime = "SURVTIME",
                          column_cnsr = "CNSR")
head(estimate_empirical_hazard(dat))
#>       ARM       time     hazard        lcl        ucl
#> 1 Control 0.07962164 0.05461771 0.04003634 0.07450965
#> 2 Control 0.50471440 0.05458969 0.04024447 0.07404828
#> 3 Control 0.54667799 0.05458692 0.04026481 0.07400339
#> 4 Control 0.54997400 0.05458671 0.04026640 0.07399987
#> 5 Control 0.72864686 0.05457494 0.04035252 0.07381011
#> 6 Control 0.73181776 0.05457473 0.04035404 0.07380676
```
